library(ncdf4)
library(PCICt)
library(proj4)

test_that("time series can be retrieved from file", {
	f <- list(t=nc_open("test1.nc"), not=nc_open("test1.nc", readunlim=FALSE))
	correct.data.ts.test1 <- structure(c(599227200, 599313600, 599400000, 599486400, 599572800),
			.Dim = 5L, cal = "365", months = c(31, 28, 31, 30, 31, 30,
					31, 31, 30, 31, 30, 31),
			class = "PCICt", dpy = 365, tzone = "GMT", units = "secs")
	correct.data.ts.test2 <- structure(c(599227200, 599313600, 599400000, 599486400, 599572800),
			.Dim = 5L, cal = "365", months = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
			class = "PCICt", dpy = 365, tzone = "GMT", units = "secs",
			bounds = structure(c(599184000, 599270400, 599270400, 599356800, 599356800, 599443200, 599443200,
							599529600, 599529600, 599616000), .Dim = c(2L, 5L), cal = "365",
					months = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), class = "PCICt", dpy = 365, tzone = "GMT", units = "secs"))
			
	lapply(f, function(f1) {
				expect_equal(nc.get.time.series(f1), correct.data.ts.test1)
				expect_equal(nc.get.time.series(f1, return.bounds=TRUE), correct.data.ts.test2)
				expect_equal(nc.get.time.series(f1, "tasmax"), correct.data.ts.test1)
				expect_equal(nc.get.time.series(f1, time.dim.name="time"), correct.data.ts.test1)
						
				expect_error(nc.get.time.series(f1, "foo"))
				expect_error(nc.get.time.series(f1, time.dim.name="foo"))
				nc_close(f1)
			})		
})

test_that("grid mappings are read and applied", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	proj4.string <- "+proj=ob_tran +o_proj=longlat +lon_0=-97 +o_lat_p=42.5 +a=1 +to_meter=0.0174532925199 +no_defs"
	expect_equal(nc.get.proj4.string(f, "tasmax"), proj4.string)
	
	lat.dat <- ncvar_get(f, "lat")
	expect_equal(dim(lat.dat), c(155,130))
	expect_equal(lat.dat[1,1], 12.3573083877563)
	expect_equal(lat.dat[1,130], 59.1459007263184)
	expect_equal(lat.dat[155,1], 12.3573160171509)
	expect_equal(lat.dat[155,130], 59.1459159851074)
	
	
	lon.dat <- ncvar_get(f, "lon")
	expect_equal(dim(lon.dat), c(155,130))
	expect_equal(lon.dat[1,1], 232.930877685547)
	expect_equal(lon.dat[1,130], 189.603073120117)
	expect_equal(lon.dat[155,1], 293.069122314453)
	expect_equal(lon.dat[155,130], 336.396881103516)
			
	indices <- matrix(c(1, 155, 155, 1, 1, 1, 130, 130), nrow=4, ncol=2)
	colnames(indices) <- c("x", "y")
	
	# The project() function in the proj4 r package returns unexpected values
	# on these inputs on Debian systems.
	# Run on Ubuntu, Windows, Fedora, CentOS, and Mac, the tests behave as expected.
	# Accordingly, the last part of this test is skipped on Debian systems.
	skip_if(grepl("Debian", Sys.info()["version"]), message="proj4 behaves unexpectedly on a Debian system")
	projected.data <- list(x=f$dim$rlon$vals[indices[,"x"]], y=f$dim$rlat$vals[indices[,"y"]])
			
	latlon.data <- project(projected.data, proj4.string, ellps.default=NA, inverse=TRUE)
	expect_equal((latlon.data$x + 360) %% 360, lon.dat[indices], tolerance=1e-5)
	expect_equal(latlon.data$y, lat.dat[indices], tolerance=1e-5)
	nc_close(f)				
})

test_that("dimensions for a compressed axis can be determined", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	compress.dims <- nc.get.compress.dims(f, "tasmax")
	expect_equal(compress.dims, list())
	nc_close(f)                                                                                                                                                                                                                                                                 
})

test_that("dimensions associated with a variable can be determined", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	dim.axes <- structure(c("Y", "X", "T", NA), .Names = c("rlat", "rlon", "time", "bnds"))
	dim.axes.var <- structure(c("X", "Y", "T"), .Names = c("rlon", "rlat", "time"))
	expect_equal(nc.get.dim.axes(f), dim.axes)
	expect_equal(nc.get.dim.axes(f, "tasmax"), dim.axes.var)
	nc_close(f)			
})

test_that("dimensions and CF axis associated with a variable can be determined", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	coord.axes <- structure(c("X", "Y"), .Names = c("lon", "lat"))
	expect_equal(nc.get.coordinate.axes(f, "tasmax"), coord.axes)
	nc_close(f)
})

test_that("CF axis of a dimension can be determined from dimension name", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	dim.axes <- structure(c(NA, NA, "T", NA), .Names = c("rlat", "rlon", "time", "bnds"))
	dim.axes.var <- structure(c(NA, NA, "T"), .Names = c("rlon", "rlat", "time"))
	expect_equal(nc.get.dim.axes.from.names(f), dim.axes)
	expect_equal(nc.get.dim.axes.from.names(f, "tasmax"), dim.axes.var)
	nc_close(f)			
})

test_that("dimension name for a file or variable can be read", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	dim.names <- c("rlat", "rlon", "time", "bnds")
	dim.names.var <- c("rlon", "rlat", "time")
	expect_equal(nc.get.dim.names(f), dim.names)
	expect_equal(nc.get.dim.names(f, "tasmax"), dim.names.var)
	nc_close(f)			
})

test_that("variables can be determined", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	expect_equal(nc.get.variable.list(f), "tasmax")
	expect_equal(nc.get.variable.list(f, min.dims=4), character(0))
	nc_close(f)
})

test_that("climatology bounds can be read", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	expect_equal(nc.get.climatology.bounds.var.list(f), logical(0))
	nc_close(f)		
})

test_that("bounds dimension can be determined", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	expect_equal(nc.get.dim.bounds.var.list(f), structure("time_bnds", .Names = "time"))
	nc_close(f)		
})

test_that("axis dimensions can be determined", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	expect_equal(nc.get.dim.for.axis(f, "tasmax", "X")$name, "rlon")
	nc_close(f)	
})

test_that("a subset of data can be read", {
	f <- nc_open("test1.nc", readunlim=FALSE)
	all.data <- ncvar_get(f, "tasmax")
	dimnames(all.data) <- list(NULL, NULL, NULL)
	expect_equal(nc.get.var.subset.by.axes(f, "tasmax", list(X=1:4, Y=c(1, 3, 5))), all.data[1:4, c(1, 3, 5), ])
	nc_close(f)
})

test_that("data can be written to a netCDF file", {
	filename <- tempfile()
	f.in <- nc_open("test1.nc")
	dat <- nc.get.var.subset.by.axes(f.in, "tasmax", list(X=1:4, Y=c(1, 3, 5)))
	x.dim <- ncdim_def("rlon", "degrees", vals=f.in$dim$rlon$vals[c(3:4, 1:2)])
	y.dim <- ncdim_def("rlat", "degrees", vals=f.in$dim$rlat$vals[c(1, 3, 5)])
	t.dim <- ncdim_def("time", "days since 1949-12-01", vals=f.in$dim$time$vals)
	var.list <- list(tasmax=ncvar_def("tasmax", "K", list(t.dim, x.dim, y.dim), 1e20, longname="Daily Maximum Near-Surface Air Temperature"))
	f.out <- nc_create(filename, var.list)
	nc.copy.atts(f.in, "rlat", f.out, "rlat")
	nc.copy.atts(f.in, "rlon", f.out, "rlon")
	nc.copy.atts(f.in, "tasmax", f.out, "tasmax")
	dat.permuted <- nc.conform.data(f.in, f.out, "tasmax", "tasmax", dat, allow.dim.subsets=TRUE)
	nc.put.var.subset.by.axes(f.out, "tasmax", dat.permuted, list())
	nc_sync(f.out)
	nc_close(f.out)
	f.out <- nc_open(filename)
			
	dat.out <- nc.get.var.subset.by.axes(f.out, "tasmax", list())
	expect_equal(as.numeric(dat.out), as.numeric(dat.permuted))
	expect_equal(dat[1, 2, 3], dat.out[3, 3, 2])
			
	nc_close(f.in)
	nc_close(f.out)
			
	unlink(filename)			
})