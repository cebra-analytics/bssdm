# Climatch SDM predict method

The model prediction component of an implementation of the ABARES
Climatch species distribution modelling (SDM) method (ABARES, 2020).

## Usage

``` r
# S3 method for class 'Climatch'
predict(
  object,
  x,
  algorithm = NULL,
  sd_data = NULL,
  as_score = NULL,
  raw_output = NULL,
  filename = "",
  parallel_cores = NULL,
  debug = FALSE,
  ...
)
```

## Arguments

- object:

  A "Climatch" model S4 object containing slots:

  `method`

  :   SDM method: "climatch".

  `algorithm`

  :   Algorithm: "euclidean" or "closest_standard_score".

  `variables`

  :   List of climate (or environmental) variable names.

  `sd`

  :   The standard deviation of each variable calculated via the climate
      data (*x*) or the *sd_data* when provided.

  `presence`

  :   The selected (nearest within range) climate data for each
      occurrence point.

  `coordinates`

  :   The coordinates for the selected climate data.

  `as_score`

  :   Indication of whether to generate a score 0-10 (TRUE) or values
      0-1 (FALSE).

- x:

  Climate (or environmental) data with corresponding model variables as
  a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  `raster::Raster*`, `data.frame`, or `matrix`.

- algorithm:

  Optional (overriding) Climatch method algorithm selected from
  "euclidean" or "closest_standard_score".

- sd_data:

  Optional (overriding) `data.frame` for calculating the standard
  deviation for climate variable, or a `vector` of pre-calculated
  values.

- as_score:

  Optional (overriding) logical to indicate whether to generate a score
  0-10 (TRUE) or values 0-1 (FALSE).

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

ABARES (2020). Climatch v2.0 User Manual. Canberra.
<https://climatch.cp1.agriculture.gov.au/> Accessed: November 2021.
