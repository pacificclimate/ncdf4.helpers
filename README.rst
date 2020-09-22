.. image:: https://github.com/pacificclimate/ncdf4.helpers/workflows/R%20CI/badge.svg
    :target: https://github.com/pacificclimate/ncdf4.helpers

.. image:: https://github.com/pacificclimate/ncdf4.helpers/workflows/R%20CI%20CRAN/badge.svg
    :target: https://github.com/pacificclimate/ncdf4.helpers


What is ncdf4.helpers?
=====================
* `ncdf4.helpers` is a collection of helpful R functions for working with netCDF files that conform to the `CF Conventions`_. The CF Conventions define metadata needed to specify the spatial and temporal properties of a dataset and a description of what each variable represents. This library can be used in addition to `ncdf4` to work with netCDF files that follow this standard.

.. _CF Conventions: https://cfconventions.org/index.html

Using ncdf4.helpers
===================

You can install ncdf4.helpers from CRAN within R ::

  > install.packages("ncdf4.helpers")
  > library("ncdf4.helpers")

Help on the package or individual functions is available within R ::

  > ?ncdf4.helpers
  > ?nc.get.time.series

 
