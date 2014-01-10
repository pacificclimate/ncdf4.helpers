ncdf4.helpers.test.file.functions <- function() {
  f1 <- nc_open("test1.nc")
  correct.data.ts.test1 <- structure(c(599227200, 599313600, 599400000, 599486400, 599572800),
                                     .Dim = 5L, cal = "365", months = c(31, 28, 31, 30, 31, 30,
                                                               31, 31, 30, 31, 30, 31),
                                     class = "PCICt", dpy = 365, tzone = "GMT", units = "secs")
  checkEquals(nc.get.time.series(f1), correct.data.ts.test1)
  
  nc_close(f1)
}
