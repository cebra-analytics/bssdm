# Range bagging SDM predict method

The model prediction component of an implementation of the range bagging
species distribution modelling (SDM) method (Drake, 2015).

## Usage

``` r
# S3 method for class 'Rangebag'
predict(
  object,
  x,
  raw_output = NULL,
  filename = "",
  parallel_cores = NULL,
  debug = FALSE,
  ...
)
```

## Arguments

- object:

  A "Rangebag" model S4 object containing slots:

  `method`

  :   SDM method: "rangebag".

  `variables`

  :   List of climate (or environmental) variable names.

  `presence`

  :   The selected climate data corresponding to occurrences points.

  `coordinates`

  :   The coordinates for the selected climate data.

  `ch_models`

  :   A list of convex hull models (vertices).

- x:

  Climate (or environmental) data with corresponding model variables as
  a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  `raster::Raster*`, `data.frame`, or `matrix`.

- raw_output:

  Logical to indicate whether to return raw predicted values (TRUE) or
  as an object (as per *x*: FALSE). Default is NULL, returning either
  raw values or a spatial raster (as per *x*).

- filename:

  Optional filename for writing spatial raster output (only). Default is
  "".

- parallel_cores:

  Optional number of cores available for parallel processing, thus
  enable parallel processing. Default is NULL (serial).

- debug:

  Output additional debug information (memory block info, etc.). Default
  is FALSE.

- ...:

  Additional parameters.

## Value

Predicted values as a raw vector or a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
`raster::Raster*`, `data.frame`, or `matrix` (as per *x*).

## References

Drake, J. M. (2015). Range bagging: a new method for ecological niche
modelling from presence-only data. *Journal of the Royal Society
Interface*, 12(107), 20150086.
[doi:10.1098/rsif.2015.0086](https://doi.org/10.1098/rsif.2015.0086)
