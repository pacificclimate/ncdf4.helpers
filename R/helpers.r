library(ncdf4)
library(PCICt)

## Returns a list of strings corresponding to bounds variables
nc.get.dim.bounds.var.list <- function(f) {
  dimension.vars <- names(f$dim)
  return(unlist(sapply(names(f$dim), function(x) {
    if(f$dim[[x]]$create_dimvar) {
      a <- ncatt_get(f, x, "bounds");
      if(a$hasatt)
        return(a$value);
    }

    ## Heuristic detection for broken files
    bnds.vars <- c(paste(x, "bnds", sep="_"), paste("bounds", x, sep="_"))
    bnds.present <- bnds.vars %in% names(f$var)
    if(any(bnds.present))
      return(bnds.vars[bnds.present])

    return(NULL);
  } )))
}

## Returns a list of strings corresponding to actual variables in files (not lat/lon/etc).
nc.get.variable.list <- function(f, con) {
  var.list <- names(f$var)
  bounds <- nc.get.dim.bounds.var.list(f)
  has.axis <- unlist(lapply(var.list, function(x) { a <- ncatt_get(f, x, "axis"); if(a$hasatt & nchar(a$value) == 1) return(x); return(NULL); } ))
  
  ## When things get really broken, we'll need this...
  bnds.heuristic <- !grepl("_bnds", var.list)
  
  var.mask <- bnds.heuristic & (!(var.list %in% c(bounds, has.axis, "lat", "lon") | unlist(lapply(f$var, function(x) { return(x$prec == "char" | x$ndims == 0) }))))
  
  return(var.list[var.mask])
}

nc.get.dim.names <- function(f, v) {
  return(unlist(lapply(f$var[[v]]$dim, function(x) { return(x$name) })))
}

nc.get.dim.axes.from.names <- function(f, v) {
  dim.names <- nc.get.dim.names(f, v)
  return(sapply(dim.names, function(x, y) { ifelse(any(x == names(y)), y[x == names(y)], NA) }, c("lat"="Y", "lon"="X", "xc"="X", "yc"="Y", "x"="X", "y"="Y", "time"="T", "plev"="Z", "lev"="Z")))
}

## Returns dimension axes according to direct (reading 'axis' attribute) and indirect (inference from dimension names) methods.
## Axes are X, Y, Z (depth, plev, etc), T (time), and S (space, for reduced grids)
nc.get.dim.axes <- function(f, v) {
  dim.names <- nc.get.dim.names(f, v)
  dim.axes <- sapply(dim.names, function(x) { if(!f$dim[[x]]$create_dimvar) return(NA); a <- ncatt_get(f, x, "axis"); return(ifelse(a$hasatt, toupper(a$value), NA)) })
  contains.compress.att <- sapply(dim.names, function(x) { ifelse(f$dim[[x]]$create_dimvar, ncatt_get(f, x, "compress")$hasatt, FALSE) })

  ## Fill in dim axes best we can if axis attributes are missing
  if(any(is.na(dim.axes)))
    dim.axes[is.na(dim.axes)] <- nc.get.dim.axes.from.names(f, v)[is.na(dim.axes)]

  if(is.list(contains.compress.att))
    browser()
  
  if(sum(contains.compress.att) != 0) 
    dim.axes[contains.compress.att] <- "S"

  return(dim.axes)
}

## Returns the dimensions used by the compressed axis
nc.get.compress.dims <- function(f, v) {
  dim.names <- nc.get.dim.names(f, v)
  dim.axes <- nc.get.dim.axes(f, v)
  compress.att <- ncatt_get(f, dim.names[dim.axes == "S"], "compress")
  compress.axes <- strsplit(compress.att$value, " ")[[1]]
  stopifnot(length(compress.axes) == 2)

  return(list(x.dim=f$dim[[which(dim.names == compress.axes[2])]], y.dim=f$dim[[which(dim.names == compress.axes[1])]]))
}

## Returns TRUE if the given data series is regular (ie: evenly spaced steps)
nc.is.regular.dimension <- function(d) {
  return(get.f.step.size(d, min) == get.f.step.size(d, max))
}

## Gets multiplier for time scale given units
nc.get.time.multiplier <- function(x) {
  return(switch(x, "days"=86400, "hours"=3600, "minutes"=60, "months"=86400 * 30))
}

## Returns the time series as PCICt
nc.get.time.series <- function(f, filename.date.parsing=FALSE, hadley.hack=FALSE, cdo.hack=FALSE, correct.for.proleptic=FALSE) {
  if(!("time" %in% names(f$dim)) || !f$dim$time$create_dimvar)
    return(NA)

  ## Hack to get around missing time axis on anomaly fields
  if(f$dim$time$len == 0) {
    if(!filename.date.parsing)
      return(NA)
    
    ## Try and infer anomaly period
    filename.bits <- strsplit(rev(strsplit(f$filename, "/")[[1]])[1], "-")[[1]]
    start.year <- filename.bits[5]
    end.year <- strsplit(filename.bits[6], "\\.")[[1]][1]
    
    return(as.PCICt(c(paste(start.year, "-01-01", sep=""), paste(end.year, "-12-31", sep="")), "gregorian"))
  }

  time.units <- f$dim$time$units
  time.split <- strsplit(f$dim$time$units, " ")[[1]]
  time.unit <- time.split[1]

  time.calendar.att <- ncatt_get(f, "time", "calendar")
  if(time.split[2] == "as") {
    ## This is to deal with retarded date formats which use format specifiers that aren't valid.
    return(as.PCICt(as.character(f$dim$time$vals), cal=ifelse(time.calendar.att$hasatt, time.calendar.att$value, "gregorian"), format=strsplit(time.split[3], "\\.")[[1]][1]))
  } else {
    time.origin.string <- time.split[3]

    if(time.split[1] == "months")
      time.origin.string <- paste(time.origin.string, "-01", sep="")
    
    if(length(time.split) > 3)
      time.origin.string <- paste(time.origin.string, time.split[4])

    cal <- ifelse(time.calendar.att$hasatt, time.calendar.att$value, "gregorian")
    cal <- ifelse(correct.for.proleptic && cal == "gregorian", "proleptic_gregorian", cal)
    
    ## Die in a fire, CDO
    ## This dirty hack is to get around CDO throwing "Gregorian" on as time units on 360-day data, but putting
    ## the time origin in the original 360-day calendar.
    if(cal != "360" & cal != "360_day" & grepl("([0-9]+)-02-30", time.origin.string))
      return(NA)

    ## FIXME: STAB HADLEY CENTER IN THE EYE
    ## Another dirty hack, this time for HadCM3:
    ## If the year field is of length 2, append "20" to it to put it in the 21st century.
    if(hadley.hack && nchar(strsplit(time.origin.string, "-")[[1]][1]) == 2)
      time.origin.string <- paste("20", time.origin.string, sep="")

    ## Specific hack for people too dumb to tell the difference between an O and a zero.
    time.origin.string <- gsub("O", "0", time.origin.string)
    
    time.origin <- as.PCICt(time.origin.string, cal=cal)
    
    time.multiplier <- nc.get.time.multiplier(time.unit)

    time.vals <- f$dim$time$vals
    if(any(is.na(time.vals)))
      time.vals <- ncvar_get(f, "time")
    
    return(time.origin + (time.vals * time.multiplier))
  }
}

get.f.step.size <- function(d, f) {
  return(match.fun(f)(d[2:length(d)] - d[1:(length(d) - 1)]))
}
