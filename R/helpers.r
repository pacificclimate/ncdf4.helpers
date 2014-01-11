library(ncdf4)
library(PCICt)
library(abind)

#' Get subsets to be distributed to workers
#'
#' Get subsets to be distributed to workers
#'
#' Given a desired number of values (\code{num.vals}), the sizes of the dimensions (\code{dim.size}), the corresponding axes (\code{dim.axes}), the desired axis to split on (\code{axis.to.split.on}), and optionally the minimum number of chunks to return (\code{min.num.chunks}), returns a list of lists of subsets appropriate to be passed to \code{nc.put.var.subsets.by.axes} or \code{nc.get.var.subsets.by.axes}.
#'
#' This functionality is useful when you want to keep memory consumption down but want to maximize the amount read in at one time to make the best use of available I/O bandwidth.
#'
#' @param num.vals The maximum number of values to process at once.
#' @param dim.size The sizes of the dimensions of the data to be processed.
#' @param dim.axes The axes of the data, as returned by \code{nc.get.dim.axes}.
#' @param axis.to.split.on The axis (X, Y, T, etc) to split the data on.
#' @param min.num.chunks The minimum number of chunks to generate, even if the chunks are considerably smaller than \code{num.vals}.
#' @return A list of lists describing subsets in a suitable form to be passed to \code{nc.put.var.subsets.by.axes} or \code{nc.get.var.subsets.by.axes}.
#'
#' @examples
#' ## Get a subset from an example
#' subsets <- get.cluster.worker.subsets(1E7, c(128, 64, 50000),
#'                                       c(lon="X", lat="Y", time="T"), "Y")
#' 
#' @export
get.cluster.worker.subsets <- function(num.vals, dim.size, dim.axes, axis.to.split.on, min.num.chunks=1) {
  split.dim <- dim.axes == axis.to.split.on
  split.dim.size <- dim.size[split.dim]
  split.slice.size <- prod(dim.size[!split.dim])
  
  rows.per.slice <- num.vals / split.slice.size
  if(split.dim.size / rows.per.slice < min.num.chunks)
    rows.per.slice <- split.dim.size / min.num.chunks
  
  stopifnot(rows.per.slice >= 1)
    
  rows.per.slice <- floor(rows.per.slice)

  row.fraction <- (0:(ceiling((split.dim.size) / rows.per.slice) - 1) * rows.per.slice)
  split.dim.starts <- floor(row.fraction) + 1
  split.dim.ends <- split.dim.starts + rows.per.slice - 1
  split.dim.ends[length(split.dim.ends)] <- split.dim.size
  return(lapply(1:length(split.dim.starts), function(x) {
    subset <- list(split.dim.starts[x]:split.dim.ends[x])
    names(subset) <- c(dim.axes[split.dim])
    subset
  }))
}

#' Splits up a CMIP5 filename.
#'
#' Splits up a CMIP5 filename into its component parts.
#'
#' This function splits up a given CMIP5 filename, returning a named vector consisting of the descriptive parts of the filename.
#'
#' @param cmip5.file The filename to be split.
#' @return A vector containing the variable (var), time resolution (tres), model (model), emissions scenario (emissions), run (run), time range (trange), time start (tstart) and time end (tend) for the file.
#'
#' @examples
#' ## Split up filename into component bits
#' split.bits <- get.split.filename.cmip5("pr/pr_day_MRI-CGCM3_historical_r1i1p1_18500101-20051231.nc")
#' 
#' @export
get.split.filename.cmip5 <- function(cmip5.file) {
  split.path <- strsplit(cmip5.file, "/")[[1]]
  fn.split <- strsplit(tail(split.path, n=1), "_")[[1]]
  names(fn.split) <- c("var", "tres", "model", "emissions", "run", "trange", rep(NA, max(0, length(fn.split) - 6)))
  fn.split[length(fn.split)] <- strsplit(fn.split[length(fn.split)], "\\.")[[1]][1]
  fn.split[c('tstart', 'tend')] <- strsplit(fn.split['trange'], "-")[[1]]
  fn.split
}

## In the case of the put method, the data is assumed to be in the same dimension order as the output file.
nc.put.subset.recursive <- function(chunked.axes.indices, f, v, dat, starts, counts, axes.map) {
  if(length(chunked.axes.indices) == 0) {
    res <- ncvar_put(f, v, dat, start=starts, count=counts)
  } else {
    axis.index <- head(chunked.axes.indices, n=1)
    axes.to.pass.on <- tail(chunked.axes.indices, n=-1)
    axis.id <- names(axis.index)

    if(length(axis.index[[1]]) == 1) {
      ## Fast path
      starts[axis.id] <- min(axis.index[[1]][[1]])
      counts[axis.id] <- length(axis.index[[1]][[1]])
      nc.put.subset.recursive(axes.to.pass.on, f, v, dat, starts, counts, axes.map)
    } else {
      ## Slow path
      cum.indices <- c(0, cumsum(sapply(axis.index, length)))
      dat.indices <- lapply(1:(length(cum.indices) - 1), function(x) { (cum.indices[x] + 1):cum.indices[x+1] })
      dat.subindex <- lapply(dim(dat), function(x) { 1:x })
      
      res <- lapply(1:length(axis.index[[1]]), function(x) {
        gc()
        indices.to.put <- axis.index[[1]][[x]]
        starts[axis.id] <- min(indices.to.put)
        counts[axis.id] <- length(indices.to.put)
        dat.subindex[[axes.map == axis.id]] <- dat.indices[[x]]
        dat.sub <- do.call('[', c(list(dat), dat.subindex, drop=FALSE))
        nc.put.subset.recursive(axes.to.pass.on, f, v, dat.sub, starts, counts, axes.map)
      })
    }
  }
}

## WARNING: Code is not fully tested.
#' Puts a data subset in the place described by the named list of axes.
#'
#' Puts a data subset in the place described by the named list of axes.
#'
#' This function will write data (\code{dat}) out to the specified file (\code{f}) and variable (\code{v}) at the location specified by \code{axis.indices}.
#'
#' @param f An object of class \code{ncdf4} which represents a NetCDF file.
#' @param v A string naming a variable in a file or an object of class \code{ncvar4}.
#' @param dat The data to put in the file.
#' @param axis.indices A list consisting of zero or more vectors of indices, named by which axis they refer to (X, Y, T, etc).
#' @param axes.map An optional vector mapping axes to NetCDF dimensions. If not supplied, it will be generated from the file.
#' @param input.axes An optional vector containing the input axis map. If supplied, it will be used to permute the data from the axis order in the input data, to the axis order in the output data.
#'
#' @examples
#' ## Copy a subset of the data from one location to another.
#' \dontrun{
#' f <- nc_open("pr.nc")
#' dat <- nc.get.var.subset.by.axes(f1, "pr", list(X=1:4, Y=c(1, 3, 5)))
#' nc.put.var.subset.by.axes(f1, "pr", dat, list(X=5:8, Y=1:3))
#' nc_close(f)
#' }
#'
#' @export
nc.put.var.subset.by.axes <- function(f, v, dat, axis.indices, axes.map=NULL, input.axes=NULL) {
  if(is.null(axes.map))
    axes.map <- nc.get.dim.axes(f, v)

  if(length(axes.map) == 0)
    return(c())

  ## Permute data to match order within file...
  if(input.axes != NULL) {
    stopifnot(length(dim(dat)) == length(input.axes))
    stopifnot(length(input.axes) == length(axes.map))
    o.axes <- order(axes.map)
    o.input <- order(input.axes)
    if(o.axes != o.input)
      dat <- aperm(dat, o.axes[o.input])
  }
  
  ## Check that all axes are in the map and that the names are the same as the dim names
  stopifnot(all(names(axis.indices) %in% axes.map))
  stopifnot(names(axes.map) %in% nc.get.dim.names(f, v))
  
  ## Chunk consecutive sets of blocks into a request
  chunked.axes.indices <- lapply(axis.indices, function(indices) {
    if(length(indices) == 0) return(c())
    boundary.indices <- c(0, which(diff(indices) != 1), length(indices))
    lapply(1:(length(boundary.indices) - 1), function(x) { (indices[boundary.indices[x] + 1]):indices[boundary.indices[x + 1]] } )
  })
  
  ## By default, fetch all data.
  starts <- rep(1, length(f$var[[v]]$dim))
  counts <- rep(-1, length(f$var[[v]]$dim))
  names(starts) <- names(counts) <- axes.map
  
  return(nc.put.subset.recursive(chunked.axes.indices, f, v, dat, starts, counts, axes.map))
}

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

    if(length(axis.index[[1]]) == 1) {
      starts[axis.id] <- min(axis.index[[1]][[1]])
      counts[axis.id] <- length(axis.index[[1]][[1]])
      nc.get.subset.recursive(axes.to.pass.on, f, v, starts, counts, axes.map)
    } else {
      res <- lapply(axis.index[[1]], function(x) {
        starts[axis.id] <- min(x)
        counts[axis.id] <- length(x)
        nc.get.subset.recursive(axes.to.pass.on, f, v, starts, counts, axes.map)
      })
      do.call(abind, c(res, list(along=which(axes.map == axis.id))))
    }
  }
}

## Why not just have a NetCDF -class- that implements a subset operator that does the actual fetching?
## WARNING: Code is not fully tested.
## FIXME: This should return the axis ordering as an attribute (or something)
## FIXME: Add a drop option to be able to replicate R's (albeit stupid) behaviour of dropping 1 length dims.
#' Gets a data subset in the place described by the named list of axes.
#'
#' Gets a data subset in the place described by the named list of axes.
#'
#' This function will read data from the specified file (\code{f}) and variable (\code{v}) at the location specified by \code{axis.indices}.
#'
#' @param f An object of class \code{ncdf4} which represents a NetCDF file.
#' @param v A string naming a variable in a file or an object of class \code{ncvar4}.
#' @param axis.indices A list consisting of zero or more vectors of indices, named by which axis they refer to (X, Y, T, etc).
#' @param axes.map An optional vector mapping axes to NetCDF dimensions. If not supplied, it will be generated from the file.
#'
#' @examples
#' ## Get a subset of the data.
#' \dontrun{
#' f <- nc_open("pr.nc")
#' dat <- nc.get.var.subset.by.axes(f1, "pr", list(X=1:4, Y=c(1, 3, 5)))
#' nc_close(f)
#' }
#'
#' @export
nc.get.var.subset.by.axes <- function(f, v, axis.indices, axes.map=NULL) {
  if(is.null(axes.map))
    axes.map <- nc.get.dim.axes(f, v)

  if(length(axes.map) == 0)
    return(c())
  
  ## Check that all axes are in the map and that the names are the same as the dim names
  stopifnot(all(names(axis.indices) %in% axes.map))
  stopifnot(names(axes.map) %in% nc.get.dim.names(f, v))
  
  ## Chunk consecutive sets of blocks into a request
  chunked.axes.indices <- lapply(axis.indices, function(indices) {
    if(length(indices) == 0) return(c())
    boundary.indices <- c(0, which(diff(indices) != 1), length(indices))
    lapply(1:(length(boundary.indices) - 1), function(x) { (indices[boundary.indices[x] + 1]):indices[boundary.indices[x + 1]] } )
  })
  
  ## By default, fetch all data.
  starts <- rep(1, length(f$var[[v]]$dim))
  counts <- rep(-1, length(f$var[[v]]$dim))
  names(starts) <- names(counts) <- axes.map
  
  return(nc.get.subset.recursive(chunked.axes.indices, f, v, starts, counts, axes.map))
}

#' Reorder data so that X and Y axes match.
#'
#' Reorder data so that X and Y axes match.
#'
#' This function will take a given current file, variable, and 3D slab of data and permute the data along the X and Y axes such that it matches the order of the data in the desired file and variable.
#'
#' @param f.desired The desired file (an object of class \code{ncdf4})
#' @param f.source The source file (an object of class \code{ncdf4})
#' @param v.desired The desired variable: a string naming a variable in a file or an object of class \code{ncvar4}.
#' @param v.source The source variable: a string naming a variable in a file or an object of class \code{ncvar4}.
#' @param dat The data to be reordered to match the XY ordering desired.
#' @return The data permuted to match the XY ordering desired.
#'
#' @examples
#' ## Copy attributes from one variable to another; but don't copy units or
#' ## standard_name, and copy long_name as old_long_name.
#' \dontrun{
#' f1 <- nc_open("pr.nc")
#' f2 <- nc_open("pr2.nc", write=TRUE)
#' dat <- nc.get.var.subset.by.axes(f1, "pr")
#' new.dat <- nc.match.xy(f2, f1, "pr", "pr", dat)
#' nc_close(f1)
#' nc_close(f2)
#' }
#'
#' @export
nc.match.xy <- function(f.desired, f.source, v.desired, v.source, dat) {
  ## Permute data to account for upside-down data and other stupid problems that cause
  ## comparisons to be irrelevant.
  ## Also subset time...
  f.desired.axes <- nc.get.dim.axes(f.desired, v.desired)
  f.source.axes <- nc.get.dim.axes(f.source, v.source)

  warning("Poorly tested function. Don't expect everything to work right.")
  
  stopifnot(f.desired.axes[1] == f.source.axes[1] && f.desired.axes[2] == f.source.axes[2])
  
  x.dim.pcic <- (f.desired$dim[[names(f.desired.axes)[f.desired.axes == "X"]]]$vals + 360) %% 360
  y.dim.pcic <- f.desired$dim[[names(f.desired.axes)[f.desired.axes == "Y"]]]$vals
  x.dim.cccma <- (f.source$dim[[names(f.source.axes)[f.source.axes == "X"]]]$vals + 360) %% 360
  y.dim.cccma <- f.source$dim[[names(f.source.axes)[f.source.axes == "Y"]]]$vals
  stopifnot(length(x.dim.pcic) == length(x.dim.cccma))
  stopifnot(length(y.dim.pcic) == length(y.dim.cccma))
  
  x.permute <- order(x.dim.cccma)[order(order(x.dim.pcic))]
  y.permute <- order(y.dim.cccma)[order(order(y.dim.pcic))]

  ## FIXME: Shaky assumptions galore here.
  dim(dat.cccma) <- dim(dat.cccma)[which(f.source.axes %in% f.desired.axes)]
  dat.cccma[x.permute, y.permute, ]
}
  
#' Copy attributes from one variable in one file to another file.
#' 
#' Copy attributes from one variable in one file to another file.
#'
#' This function copies attributes from a variable in one file to a variable in another file. If the source or destination variable is 0, then attributes are copied from/to the NetCDF file's global attributes.
#'
#' If desired, some attributes can be left out using \code{exception.list}, a vector of names of attributes to be excluded.
#' 
#' Attributes can also be renamed at the destination using \code{rename.mapping}, a named vector of strings in which the name of the attribute to be renamed is the name, and the attribute's new name is the value.
#'
#' @param f.src The source file (an object of class \code{ncdf4})
#' @param v.src The source variable: a string naming a variable in a file or an object of class \code{ncvar4}.
#' @param f.dest The destination file (an object of class \code{ncdf4})
#' @param v.dest The destination variable: a string naming a variable in a file or an object of class \code{ncvar4}.
#' @param exception.list A vector containing names of variables not to be copied.
#' @param rename.mapping A vector containing named values mapping source to destination names.
#' @param definemode Whether the file is already in define mode.
#'
#' @examples
#' ## Copy attributes from one variable to another; but don't copy units or
#' ## standard_name, and copy long_name as old_long_name.
#' \dontrun{
#' f1 <- nc_open("pr.nc")
#' f2 <- nc_open("pr2.nc")
#' nc.copy.atts(f1, "pr", f2, "pr", c("units", "standard_name"),
#'              c(long_name="old_long_name"))
#' dim.axes <- nc.get.dim.axes.from.names(f, "pr")
#' nc_close(f1)
#' nc_close(f2)
#' }
#'
#' @export
nc.copy.atts <- function(f.src, v.src, f.dest, v.dest, exception.list=NULL, rename.mapping=NULL, definemode=FALSE) {
  atts <- ncatt_get(f.src, v.src)
  if(length(atts) > 0) {
    if(!definemode)
      nc_redef(f.dest)

    lapply(names(atts), function(x) {
      if(!(x %in% exception.list)) {
        att.name <- if(x %in% names(rename.mapping)) rename.mapping[x] else x
        ncatt_put(f.dest, v.dest, att.name, atts[[x]], definemode=TRUE)
      }
    })

    if(!definemode)
      nc_enddef(f.dest)
  }
}

## Returns the values corresponding to the dimension variable in question
#' Get dimension corresponding to a given axis
#'
#' Get dimension corresponding to a given axis
#'
#' This function returns the dimension (of class 'ncdim4') corresponding to the specified axis (X, Y, Z, T, or S).
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v The source variable: a string naming a variable in a file or an object of class \code{ncvar4}.
#' @param axis The axis to retrieve the dimension for: a string consisting of either X, Y, Z, T, or S.
#' @return An object of class \code{ncdim4} if a dimension is found for the specified axis; \code{NA} otherwise.
#'
#' @examples
#' ## Get dimension for X axis
#' \dontrun{
#' f <- nc_open("pr.nc")
#' x.axis.dim <- nc.get.dim.axes.from.names(f, "pr", "X")
#' nc_close(f)
#' }
#'
#' @export
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
#' Get a list of names of dimension bounds variables.
#' 
#' Get a list of names of dimension bounds variables.
#'
#' This function returns the names of any dimension bounds variables found in a file.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @return A character vector naming all of the dimension bounds variables found.
#'
#' @examples
#' ## Get list of dimension bound variables
#' \dontrun{
#' f <- nc_open("pr.nc")
#' dim.bounds.var.list <- nc.get.dim.bounds.var.list(f)
#' nc_close(f)
#' }
#'
#' @export
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
#' Get a list of names of climatology bounds variables.
#' 
#' Get a list of names of climatology bounds variables.
#'
#' This function returns the names of any climatology bounds variables found in a file.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @return A character vector naming all of the climatology bounds variables found.
#'
#' @examples
#' ## Get list of climatology bounds variables
#' \dontrun{
#' f <- nc_open("pr.nc")
#' dim.axes <- nc.get.climatology.bounds.var.list(f)
#' nc_close(f)
#' }
#'
#' @export
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
#' Get a list of names of data variables.
#' 
#' Get a list of names of data variables.
#'
#' This function returns the names of any data variables found in the file -- that is, variables which are NOT dimension variables, dimension bounds variables, climatology bounds variables, coordinate variables, or grid mapping variables.
#' 
#' Optionally, one may require that the variables have a minimum number of dimensions; this can eliminate unwanted variables left in files.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param min.dims The minimum number of dimensions a variable must have to be included.
#' @return A character vector naming all of the data variables found.
#'
#' @examples
#' ## Get dimension axes from file by inferring them from dimension names
#' \dontrun{
#' f <- nc_open("pr.nc")
#' var.list <- nc.get.variable.list(f)
#' nc_close(f)
#' }
#'
#' @export
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

#' Get a list of names of dimensions.
#' 
#' Get a list of names of dimensions.
#'
#' This function returns the names of dimensions in a file or, if \code{v} is also supplied, attached to a particular variable.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v Optionally, a variable
#' @return A character vector naming the dimensions found.
#'
#' @examples
#' ## Get dimension names
#' \dontrun{
#' f <- nc_open("pr.nc")
#' dim.names <- nc.get.dim.names(f, "pr")
#' nc_close(f)
#' }
#'
#' @export
nc.get.dim.names <- function(f, v) {
  if(missing(v)) {
    d <- unlist(lapply(f$dim, function(x) { return(x$name) }))
    names(d) <- NULL
    return(d)
  } else
    return(unlist(lapply(f$var[[v]]$dim, function(x) { return(x$name) })))
}

#' Infer dimension axes from names of dimensions.
#' 
#' Infer dimension axes from names of dimensions.
#'
#' This function makes educated guesses as to what axes dimensions may apply to in the case of files with poor metadata.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v The name of a variable
#' @param dim.names Optionally, dimension names (to avoid looking them up repeatedly)
#' @return A named character vector mapping dimension names to axes.
#'
#' @examples
#' ## Get dimension axes from file by inferring them from dimension names
#' \dontrun{
#' f <- nc_open("pr.nc")
#' dim.axes <- nc.get.dim.axes.from.names(f, "pr")
#' nc_close(f)
#' }
#'
#' @export
nc.get.dim.axes.from.names <- function(f, v, dim.names) {
  if(missing(dim.names))
    dim.names <- nc.get.dim.names(f, v)
  return(sapply(dim.names, function(x, y) { ifelse(any(x == names(y)), y[x == names(y)], NA) }, c("lat"="Y", "latitude"="Y", "lon"="X", "longitude"="X", "xc"="X", "yc"="Y", "x"="X", "y"="Y", "time"="T", "timeofyear"="T", "plev"="Z", "lev"="Z", "level"="Z")))
}

#' Get a list of dimension variables and axes for a variable's coordinate variable
#' 
#' Get a list of dimension variables and axes for a variable's coordinate variable
#'
#' This function returns a named list of axes, the names of which are the associated dimension variables.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v The name of a variable
#' @return A named character vector containing axes, the names of which are the corresponding dimension variables.
#'
#' @examples
#' ## Get coordinate axes from file.
#' \dontrun{
#' f <- nc_open("pr.nc")
#' coord.axes <- nc.get.coordinate.axes(f, "pr")
#' nc_close(f)
#' }
#'
#' @export
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
#' Get dimension axes
#' 
#' Get dimension axes for the given variable
#'
#' This function returns the dimension axes for a given variable as a named character vector; the names are the names of the corresponding dimensions. If no variable is supplied, the function will return data for all dimensions found in the file.
#'
#' Axes are X, Y, Z (depth, plev, etc), T (time), and S (space, for reduced grids).
#'
#' This routine will attempt to infer axes for dimensions if no 'axis' attribute is found on a dimension variable, using the nc.get.dim.axes.from.names function.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v The name of a variable
#' @param dim.names Optionally, dimension names (to avoid looking them up repeatedly)
#' @return A named character vector mapping dimension names to axes.
#'
#' @examples
#' ## Get dimension axes from file.
#' \dontrun{
#' f <- nc_open("pr.nc")
#' ## Get dim axes for a specified variable
#' dim.axes <- nc.get.dim.axes(f, "pr")
#' ## Get all dim axes in file
#' dim.axes <- nc.get.dim.axes(f)
#' nc_close(f)
#' }
#'
#' @export
nc.get.dim.axes <- function(f, v, dim.names) {
  if(missing(dim.names))
    if(missing(v))
      dim.names <- nc.get.dim.names(f)
    else
      dim.names <- nc.get.dim.names(f, v)

  if(length(dim.names) == 0)
    return(c())

  has.dim.no.data <- function(x) { !is.null(f$dim[[x]]) && !is.null(f$dim[[x]]$create_dimvar) && !f$dim[[x]]$create_dimvar }
  
  dim.axes <- sapply(dim.names, function(x) { if(has.dim.no.data(x)) return(NA); a <- ncatt_get(f, x, "axis"); return(ifelse(a$hasatt, toupper(a$value), NA)) })
  contains.compress.att <- sapply(dim.names, function(x) { ifelse(has.dim.no.data(x) || is.null(f$var[[x]]), FALSE, ncatt_get(f, x, "compress")$hasatt) })

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
#' Get X and Y dimension variables for reduced (compressed) grids.
#' 
#' Get X and Y dimension variables for reduced (compressed) grids.
#'
#' This function retrieves the X and Y dimensions for reduced (compressed) grids, returning a list containing the X and Y dimensions.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v The name of a variable
#' @return A list consisting of two members of class \code{ncdim4}: x.dim for the X axis, and y.dim for the Y axis.
#'
#' @examples
#' ## Get compress dimensions from file.
#' \dontrun{
#' f <- nc_open("pr.nc")
#' compress.dims <- nc.get.compress.dims(f, "pr")
#' nc_close(f)
#' }
#'
#' @export
nc.get.compress.dims <- function(f, v) {
  dim.names <- nc.get.dim.names(f, v)
  dim.axes <- nc.get.dim.axes(f, v)
  if(sum(dim.axes == "S", na.rm=TRUE) == 0)
    return(list())
  compress.att <- ncatt_get(f, dim.names[dim.axes == "S"], "compress")
  compress.axes <- strsplit(compress.att$value, " ")[[1]]
  stopifnot(length(compress.axes) == 2)

  return(list(x.dim=f$dim[[which(dim.names == compress.axes[2])]], y.dim=f$dim[[which(dim.names == compress.axes[1])]]))
}

## Returns TRUE if the given data series is regular (ie: evenly spaced steps)
## Tolerance is as a fraction
#' Determine if a dimension is regular
#' 
#' Determine if a dimension is regular
#' 
#' Given supplied data and optionally a tolerance level, determine if the dimension is regular or not.
#'
#' @param d The data to be tested
#' @param tolerance The tolerance for variation in step size.
#' @return TRUE if the data is regular; FALSE if not.
#'
#' @examples
#' dat <- c(1, 2, 3, 4, 5, 6, 7)
#' ## TRUE
#' nc.is.regular.dimension(dat)
#'
#' dat[7] <- 7.001
#' ## FALSE
#' nc.is.regular.dimension(dat)
#' 
#' @export
nc.is.regular.dimension <- function(d, tolerance=0.000001) {
  return(abs((get.f.step.size(d, min) / get.f.step.size(d, max)) - 1) < tolerance)
}

#' Gets multiplier for time scale given units
#' 
#' Gets multiplier for time scale given units
#'
#' Given supplied units (days, hours, etc) returns a multiplier to convert the units into seconds.
#'
#' @param x The time scale
#' @return A numeric conversion factor to convert to seconds.
#'
#' @examples
#' ## Will return 3600
#' mul <- nc.get.time.multiplier("hours")
#'
#' @export
nc.get.time.multiplier <- function(x) {
  return(switch(x, "days"=86400, "hours"=3600, "minutes"=60, "months"=86400 * 30))
}

## Returns the time series as PCICt
#' Returns time axis data as PCICt for a file
#'
#' Returns time axis data as PCICt for a file
#'
#' This function returns time data for a file as PCICt, doing all necessary conversions.
#'
#' @note If the file was opened with \code{readunlim=FALSE}, it will read in the time values from the file; otherwise, it will retrieve the time values from the \code{ncdf4} class' data structures.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v Optionally, the variable to look for a time dimension on.
#' @param time.dim.name Optionally, the time dimension name.
#' @param correct.for.gregorian.julian Specific workaround for Gregorian-Julian calendar transitions in non-proleptic Gregorian calendars
#' @param return.bounds Whether to return the time bounds as an additional attribute
#' @return A vector of PCICt objects, optionally with bounds
#'
#' @examples
#' ## Get time series from file
#' \dontrun{
#' f <- nc_open("pr.nc")
#' ts <- nc.get.time.series(f)
#' nc_close(f)
#' }
#'
#' @export
nc.get.time.series <- function(f, v, time.dim.name, correct.for.gregorian.julian=FALSE, return.bounds=FALSE) {
  ## If the time dim wasn't supplied, go find it.
  if(missing(time.dim.name)) {
    if(missing(v))
      dim.axes <- nc.get.dim.axes(f)
    else {
      if(!(v %in% names(f$var)))
        stop(paste("Variable '", v, "' not found in file.", sep=""))
      dim.axes <- nc.get.dim.axes(f, v)
    }
    
    num.T.axes <- sum(dim.axes == "T", na.rm=TRUE)
    if(num.T.axes == 0)
      return(NA)
    else if(num.T.axes > 1)
      stop("More than one time axis found; please specify a variable or provide a name for the time axis.")

    time.dim.name <- names(dim.axes[dim.axes == "T" & !is.na(dim.axes)])
  }

  if(!(time.dim.name %in% names(f$dim)))
    stop(paste("Couldn't find dimension '", time.dim.name, "' in file.", sep=""))

  if(!f$dim[[time.dim.name]]$create_dimvar)
    stop(paste("Couldn't find dimension variable for dim '", time.dim.name, "' in file.", sep=""))

  if(f$dim$time$len == 0) {
    return(NA)
  }

  time.units <- f$dim$time$units
  time.split <- strsplit(f$dim$time$units, " ")[[1]]
  time.res <- time.split[1]

  time.calendar.att <- ncatt_get(f, time.dim.name, "calendar")
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
    
    ## Specific hack for people too dumb to tell the difference between an O and a zero.
    time.origin.string <- gsub("O", "0", time.origin.string)
    
    time.origin <- as.PCICt(time.origin.string, cal=cal)
    
    time.multiplier <- nc.get.time.multiplier(time.res)

    time.vals <- f$dim$time$vals
    if(any(is.na(time.vals)))
      time.vals <- ncvar_get(f, time.dim.name)

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

    ## Bounds processing
    bounds.vals <- NULL
    if(return.bounds) {
      bounds.att <- ncatt_get(f, time.dim.name, "bounds")
      if(bounds.att$hasatt) {
        bounds.vals <- ncvar_get(f, bounds.att$value)
      }
    }

    return(if(return.bounds) structure(time.origin + julian.correction + (time.vals * time.multiplier), bounds=time.origin + (bounds.vals * time.multiplier)) else structure(time.origin + julian.correction + (time.vals * time.multiplier)))
  }
}

#' Creates time bounds for a time series
#'
#' Creates time bounds for a time series
#'
#' Given a time series of PCICt, returns a set of bounds for that time series based the supplied units.
#'
#' @param ts The time values, of type \code{PCICt}
#' @param unit The units to be used.
#' @return 2-dimensional bounds array for the time values with dimensions [length(ts), 2].
#'
#' @examples
#' ts <- as.PCICt(c("1961-01-15", "1961-02-15", "1961-03-15"), cal="360")
#' ts.bounds <- nc.make.time.bounds(ts, unit="month")
#'
#' @export
nc.make.time.bounds <- function(ts, unit=c("year", "month")) {
  unit <- match.arg(unit)
  multiplier <- switch(unit, year=1, month=12)
  r <- range(ts)
  r.years <- as.numeric(format(r, "%Y"))
  start.date <- as.PCICt(paste(r.years[1], "-01-01", sep=""), attr(ts, "cal"))
  num.years <- r.years[2] - r.years[1] + 1
  padded.dates <- seq(start.date, by=paste("1", unit), length.out=num.years * multiplier + 1)
  padded.length <- length(padded.dates)
  bounds <- c(padded.dates[1:(padded.length - 1)], padded.dates[2:padded.length] - 86400)
  dim(bounds) <- c(padded.length - 1, 2)
  t(bounds)
}

##nc.apply <- function(var, margin, fun, nc.file, chunk.size.mb=1000, ...) {
##  if(any(margin > length(nc.file$var[[var]]$dim
##  
##}

#' Get step size for data
#'
#' Get step size for data
#'
#' Gets the step size for data, aggregated by the supplied function.
#'
#' @param d The data to have the step size determined
#' @param f The function to aggregate the step size
#' @return The step size
#'
#' @examples
#' dat <- c(1, 2, 3, 4, 5, 7)
#' max.step.size <- get.f.step.size(dat, max)
#' min.step.size <- get.f.step.size(dat, min)
#'
#' @export
get.f.step.size <- function(d, f) {
  return(match.fun(f)(diff(d, lag=1)))
}

normalize180 <- function(x) {
  (x + 180) %% 360 - 180
}

nc.get.polar.stereo.proj4.string <- function(f, grid.mapping.name) {
  lat.ts.att <- ncatt_get(f, grid.mapping.name, "standard_parallel")
  lat.0.att <- ncatt_get(f, grid.mapping.name, "latitude_of_projection_origin")
  lon.0.att <- ncatt_get(f, grid.mapping.name, "straight_vertical_longitude_from_pole")
  x.0.att <- ncatt_get(f, grid.mapping.name, "false_easting")
  y.0.att <- ncatt_get(f, grid.mapping.name, "false_northing")

  stopifnot(lat.ts.att$hasatt & lat.0.att$hasatt & lon.0.att$hasatt & x.0.att$hasatt & y.0.att$hasatt)

  if(x.0.att$value == "")
    x.0.att$value <- 0
  
  if(y.0.att$value == "")
    y.0.att$value <- 0
  
  return(paste("+proj=stere +lat_ts=", lat.ts.att$value, " +lat_0=", lat.0.att$value, " +lon_0=", lon.0.att$value, " +x_0=", x.0.att$value, " +y_0=", y.0.att$value, " +k_0=1", sep=""))
}

nc.get.rotated.pole.proj4.string <- function(f, grid.mapping.name) {
  lat.0.att <- ncatt_get(f, grid.mapping.name, "north_pole_latitude")
  lon.0.att <- ncatt_get(f, grid.mapping.name, "north_pole_longitude")
  if(!(lat.0.att$hasatt & lon.0.att$hasatt)) {
    lat.0.att <- ncatt_get(f, grid.mapping.name, "grid_north_pole_latitude")
    lon.0.att <- ncatt_get(f, grid.mapping.name, "grid_north_pole_longitude")
  }
  stopifnot(lat.0.att$hasatt & lon.0.att$hasatt)

  ## The more or less direct way here is to generate an inverse projection by feeding the values directly in as o_lon_p and o_lat_p; this is to generate a normal, forward projection.
  return(paste("+proj=ob_tran +o_proj=latlon +lon_0=", normalize180(lon.0.att$value + 180), " +o_lat_p=", lat.0.att$value, " +a=1 +to_meter=0.0174532925199 +no_defs", sep=""))
}

nc.get.lambert.conformal.conic.proj4.string <- function(f, grid.mapping.name) {
  lat.ts.att <- ncatt_get(f, grid.mapping.name, "standard_parallel")

  lat.0.att <- ncatt_get(f, grid.mapping.name, "latitude_of_projection_origin")
  lon.0.att <- ncatt_get(f, grid.mapping.name, "longitude_of_central_meridian")
  x.0.att <- ncatt_get(f, grid.mapping.name, "false_easting")
  y.0.att <- ncatt_get(f, grid.mapping.name, "false_northing")

  stopifnot(lat.ts.att$hasatt & lat.0.att$hasatt & lon.0.att$hasatt & x.0.att$hasatt & y.0.att$hasatt)
  
  return(paste("+proj=lcc +lat_0=", lat.0.att$value, " +lat_1=", lat.ts.att$value[1], " +lat_2=", lat.ts.att$value[2], " +lon_0=", lon.0.att$value, " +y_0=", y.0.att$value, " +x_0=", x.0.att$value, sep=""))
}

nc.get.transverse.mercator.proj4.string <- function(f, grid.mapping.name) {
  lat.0.att <- ncatt_get(f, grid.mapping.name, "latitude_of_projection_origin")
  lon.0.att <- ncatt_get(f, grid.mapping.name, "longitude_of_central_meridian")
  k.0.att <- ncatt_get(f, grid.mapping.name, "scale_factor_at_central_meridian")
  x.0.att <- ncatt_get(f, grid.mapping.name, "false_easting")
  y.0.att <- ncatt_get(f, grid.mapping.name, "false_northing")

  stopifnot(k.0.att$hasatt & lat.0.att$hasatt & lon.0.att$hasatt & x.0.att$hasatt & y.0.att$hasatt)

  return(paste("+proj=tmerc +lat_0=", lat.0.att$value, " +lon_0=", lon.0.att$value, " +k_0=", k.0.att$value, " +y_0=", y.0.att$value, " +x_0=", x.0.att$value, sep=""))
}

## Returns the spatial reference ID of the data set, or WGS84 (4326) if nothing found
#' Gets the proj4 string for a file
#'
#' Gets the proj4 string for a file
#'
#' Given a file and a variable, attempts to determine what the proj4 string for the given file should be. If no projection data is found, returns an empty string. Currently supports Lambert Conformal Conic, Transverse Mercator, Polar Sterographic, and Rotated Pole projections.
#'
#' @param f The file (an object of class \code{ncdf4})
#' @param v The name of a variable
#' @return A string containing the proj4 string
#'
#' @examples
#' ## Get the proj4 string for a hypothetical file.
#' \dontrun{
#' f <- nc_open("pr.nc")
#' proj4.string <- nc.get.proj4.string(f, "pr")
#' nc_close(f)
#' }
#'
#' @export
nc.get.proj4.string <- function(f, v) {
  grid.mapping.att <- ncatt_get(f, v, "grid_mapping")
  if(!grid.mapping.att$hasatt) {
    return("");
  } else {
    grid.mapping.name.att <- ncatt_get(f, grid.mapping.att$value, "grid_mapping_name")
    
    proj4.string <- switch(grid.mapping.name.att$value,
                           polar_stereographic=nc.get.polar.stereo.proj4.string(f, grid.mapping.att$value),
                           rotated_latitude_longitude=nc.get.rotated.pole.proj4.string(f, grid.mapping.att$value),
                           lambert_conformal_conic=nc.get.lambert.conformal.conic.proj4.string(f, grid.mapping.att$value),
                           transverse_mercator=nc.get.transverse.mercator.proj4.string(f, grid.mapping.att$value)
                           )
    return(proj4.string)
  }
}

