library(ncdf4)
library(PCICt)
library(abind)

nc.get.subset.recursive <- function(chunked.axes.indices, f, v, starts, counts, axes.map) {
  if(length(chunked.axes.indices) == 0) {
    res <- ncvar_get(f, v, start=starts, count=counts, collapse_degen=FALSE)
    ## Work around dwpierce's dropping of 1-dim dims in input.
    res.dim <- counts
    res.dim[counts == -1] <- f$var[[v]]$varsize[counts == -1]
    dim(res) <- res.dim
    res
  } else {
    axis.index <- head(chunked.axes.indices, n=1)
    axes.to.pass.on <- tail(chunked.axes.indices, n=-1)
    axis.id <- names(axis.index)

    res <- lapply(axis.index[[1]], function(x) {
      starts[axis.id] <- min(x)
      counts[axis.id] <- length(x)
      nc.get.subset.recursive(axes.to.pass.on, f, v, starts, counts, axes.map)
    })
    do.call(abind, c(res, list(along=which(axes.map == axis.id))))
  }
}

## Why not just have a NetCDF -class- that implements a subset operator that does the actual fetching?
## WARNING: Code is not fully tested.
nc.get.var.subset.by.axes <- function(f, v, axis.indices, axes.map=NULL) {
  if(is.null(axes.map))
    axes.map <- nc.get.dim.axes(f, v)
  
  ## Check that all axes are in the map and that the names are the same as the dim names
  stopifnot(all(names(axis.indices) %in% axes.map))
  stopifnot(names(axes.map) %in% nc.get.dim.names(f, v))
  
  ## Chunk consecutive sets of blocks into a request
  chunked.axes.indices <- lapply(axis.indices, function(indices) {
    if(length(indices) == 0) return(c())
    boundary.indices <- c(0, which(diff(indices) != 1), length(indices))
    boundary.matrix <- rbind(boundary.indices[1:(length(boundary.indices) - 1)] + 1, boundary.indices[2:length(boundary.indices)])
    as.data.frame(apply(boundary.matrix, 2, function(x) { indices[x[1]:x[2]] } ))
  })
  
  ## By default, fetch all data.
  starts <- rep(1, length(f$var[[v]]$dim))
  counts <- rep(-1, length(f$var[[v]]$dim))
  names(starts) <- names(counts) <- axes.map
  
  return(nc.get.subset.recursive(chunked.axes.indices, f, v, starts, counts, axes.map))
}

nc.permute.data.to.match <- function(f.pcic, f.cccma, v.pcic, v.cccma, dat.cccma) {
  ## Permute data to account for upside-down data and other stupid problems that cause
  ## comparisons to be irrelevant.
  ## Also subset time...
  f.pcic.axes <- nc.get.dim.axes(f.pcic, v.pcic)
  f.cccma.axes <- nc.get.dim.axes(f.cccma, v.cccma)

  warning("Poorly tested function. Don't expect everything to work right.")
  
  stopifnot(f.pcic.axes[1] == f.cccma.axes[1] && f.pcic.axes[2] == f.cccma.axes[2])
  
  x.dim.pcic <- (f.pcic$dim[[names(f.pcic.axes)[f.pcic.axes == "X"]]]$vals + 360) %% 360
  y.dim.pcic <- f.pcic$dim[[names(f.pcic.axes)[f.pcic.axes == "Y"]]]$vals
  x.dim.cccma <- (f.cccma$dim[[names(f.cccma.axes)[f.cccma.axes == "X"]]]$vals + 360) %% 360
  y.dim.cccma <- f.cccma$dim[[names(f.cccma.axes)[f.cccma.axes == "Y"]]]$vals
  stopifnot(length(x.dim.pcic) == length(x.dim.cccma))
  stopifnot(length(y.dim.pcic) == length(y.dim.cccma))
  
  x.permute <- order(x.dim.cccma)[order(order(x.dim.pcic))]
  y.permute <- order(y.dim.cccma)[order(order(y.dim.pcic))]

  ## FIXME: Shaky assumptions galore here.
  dim(dat.cccma) <- dim(dat.cccma)[which(f.cccma.axes %in% f.pcic.axes)]
  dat.cccma[x.permute, y.permute, ]
}
  
## Copy attributes from one variable in one file, to another file
nc.copy.atts <- function(f.src, v.src, f.dest, v.dest, exception.list=NULL, definemode=FALSE) {
  atts <- ncatt_get(f.src, v.src)
  if(length(atts) > 0) {
    if(!definemode)
      nc_redef(f.dest)

    ## Copy atts, with or without an exception list
    if(is.null(exception.list)) {
      lapply(names(atts), function(x) {
        ncatt_put(f.dest, v.dest, x, atts[[x]], definemode=TRUE)
      })
    } else {
      lapply(names(atts), function(x) {
        if(!(x %in% exception.list))
          ncatt_put(f.dest, v.dest, x, atts[[x]], definemode=TRUE)
      })
    }

    if(!definemode)
      nc_enddef(f.dest)
  }
}

## Returns the values corresponding to the dimension variable in question
nc.get.dim.for.axis <- function(f, v, axis) {
  dims <- f$var[[v]]$dim
  axes <- nc.get.dim.axes(f, v)

  axis.number <- which(axes == axis)
  if(length(axis.number) == 1) {
    return(dims[[axis.number]])
  } else {
    return(NA)
  }
}

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

## Returns a list of climatology bounds variables.
nc.get.climatology.bounds.var.list <- function(f) {
  dim.list <- names(f$dim)
  is.climatology<- sapply(dim.list, function(x) {
    if(f$dim[[x]]$create_dimvar && f$dim[[x]]$unlim) {
      a <- ncatt_get(f, x, "climatology")
      if(a$hasatt)
        return(a$value)
    }
    return(NA)
  })
  return(unique(is.climatology[!is.na(is.climatology)]))
}

## Returns a list of strings corresponding to actual variables in files (not lat/lon/etc).
nc.get.variable.list <- function(f, min.dims=1) {
  var.list <- names(f$var)
  enough.dims <- sapply(var.list, function(v) { length(f$var[[v]]$dim) >= min.dims } )
  bounds <- nc.get.dim.bounds.var.list(f)
  climatology.bounds <- nc.get.climatology.bounds.var.list(f)
  has.axis <- unlist(lapply(var.list, function(x) { a <- ncatt_get(f, x, "axis"); if(a$hasatt & nchar(a$value) == 1) return(x); return(NULL); } ))
  
  ## When things get really broken, we'll need this...
  bnds.heuristic <- !grepl("_bnds", var.list)
  
  var.mask <- bnds.heuristic & enough.dims & (!(var.list %in% c(bounds, has.axis, climatology.bounds, "lat", "lon") | unlist(lapply(f$var, function(x) { return(x$prec == "char" | x$ndims == 0) }))))
  
  return(var.list[var.mask])
}

nc.get.dim.names <- function(f, v) {
  return(unlist(lapply(f$var[[v]]$dim, function(x) { return(x$name) })))
}

nc.get.dim.axes.from.names <- function(f, v, dim.names) {
  if(missing(dim.names))
    dim.names <- nc.get.dim.names(f, v)
  return(sapply(dim.names, function(x, y) { ifelse(any(x == names(y)), y[x == names(y)], NA) }, c("lat"="Y", "latitude"="Y", "lon"="X", "longitude"="X", "xc"="X", "yc"="Y", "x"="X", "y"="Y", "time"="T", "timeofyear"="T", "plev"="Z", "lev"="Z", "level"="Z")))
}

nc.get.coordinate.axes <- function(f, v) {
  coords.att <- ncatt_get(f, v, "coordinates")
  if(coords.att$hasatt) {
    split.bits <- strsplit(coords.att$value, " ")[[1]]
    coords.axes <- nc.get.dim.axes(f, dim.names=split.bits)
    return(coords.axes)
  } else {
    return(c())
  }
}

## Returns dimension axes according to direct (reading 'axis' attribute) and indirect (inference from dimension names) methods.
## Axes are X, Y, Z (depth, plev, etc), T (time), and S (space, for reduced grids)
nc.get.dim.axes <- function(f, v, dim.names) {
  if(missing(dim.names))
     dim.names <- nc.get.dim.names(f, v)
     
  dim.axes <- sapply(dim.names, function(x) { if((!is.null(f$dim[[x]]) && !f$dim[[x]]$create_dimvar) || is.null(f$var[[x]])) return(NA); a <- ncatt_get(f, x, "axis"); return(ifelse(a$hasatt, toupper(a$value), NA)) })
  contains.compress.att <- sapply(dim.names, function(x) { ifelse((!is.null(f$dim[[x]]) && f$dim[[x]]$create_dimvar) || !is.null(f$var[[x]]), ncatt_get(f, x, "compress")$hasatt, FALSE) })

  ## Fill in dim axes best we can if axis attributes are missing
  if(any(is.na(dim.axes)))
    dim.axes[is.na(dim.axes)] <- nc.get.dim.axes.from.names(f, v, dim.names)[is.na(dim.axes)]

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
nc.get.time.series <- function(f, filename.date.parsing=FALSE, hadley.hack=FALSE, cdo.hack=FALSE, correct.for.gregorian.julian=FALSE, return.bounds=FALSE) {
  ## FIXME: Identify dim by axis here...
  time.var.name <- "time"
  if(!(time.var.name %in% names(f$dim)) || !f$dim$time$create_dimvar)
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

  time.calendar.att <- ncatt_get(f, time.var.name, "calendar")
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
    
    ## Die in a fire, CDO
    ## This dirty hack is to get around CDO throwing "Gregorian" on as time units on 360-day data, but putting
    ## the time origin in the original 360-day calendar.
    if(cdo.hack & cal != "360" & cal != "360_day" & grepl("([0-9]+)-02-30", time.origin.string))
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

    ## Bounds processing
    bounds.vals <- NULL
    if(return.bounds) {
      bounds.att <- ncatt_get(f, time.var.name, "bounds")
      if(bounds.att$hasatt) {
        bounds.vals <- ncvar_get(f, bounds.att$value)
      }
    }
        
    time.vals <- f$dim$time$vals
    if(any(is.na(time.vals)))
      time.vals <- ncvar_get(f, time.var.name)

    ## Correct calendar output if true gregorian
    seconds.per.day <- 86400
    origin.year.POSIXlt <- 1900
    x <- as.POSIXlt(time.origin)
    julian.correction <- 0
    if(correct.for.gregorian.julian && cal == "gregorian") {
      year.adjusted <- x$year + origin.year.POSIXlt + as.numeric(x$mon >= 3) - 1
      diff.days <- floor(year.adjusted / 100) - floor(year.adjusted / 400) - 2
      diff.days[x$year > 1582 | (x$year == 1582 & (x$mon > 10 | (x$mon == 10 & x$day >= 4))) ] <- 0
      julian.correction <- diff.days * seconds.per.day
    }
    
    return(structure(time.origin + julian.correction + (time.vals * time.multiplier), bounds=time.origin + (bounds.vals * time.multiplier)))
  }
}

##nc.apply <- function(var, margin, fun, nc.file, chunk.size.mb=1000, ...) {
##  if(any(margin > length(nc.file$var[[var]]$dim
##  
##}

get.f.step.size <- function(d, f) {
  return(match.fun(f)(d[2:length(d)] - d[1:(length(d) - 1)]))
}
