# Climatch SDM model building method

The model building component of an implementation of the ABARES Climatch
species distribution modelling (SDM) method (ABARES, 2020).

## Usage

``` r
climatch(
  x,
  p,
  algorithm = c("euclidean", "closest_standard_score"),
  d_max = 50,
  sd_data = NULL,
  as_score = TRUE,
  ...
)

# S3 method for class 'Raster'
climatch(
  x,
  p,
  algorithm = c("euclidean", "closest_standard_score"),
  d_max = 50,
  sd_data = NULL,
  as_score = TRUE,
  ...
)

# S3 method for class 'SpatRaster'
climatch(
  x,
  p,
  algorithm = c("euclidean", "closest_standard_score"),
  d_max = 50,
  sd_data = NULL,
  as_score = TRUE,
  ...
)

# S3 method for class 'SpatVector'
climatch(
  x,
  p,
  algorithm = c("euclidean", "closest_standard_score"),
  d_max = 50,
  sd_data = NULL,
  as_score = TRUE,
  ...
)

# S3 method for class 'data.frame'
climatch(
  x,
  p,
  algorithm = c("euclidean", "closest_standard_score"),
  d_max = 50,
  sd_data = NULL,
  as_score = TRUE,
  ...
)
```

## Arguments

- x:

  Climate (or environmental) data as a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  or `raster::Raster*` (with any CRS), or a `data.frame` with WGS84
  *lon* and *lat* columns.

- p:

  Species occurrence data as a `data.frame` (or `matrix`) with WGS84
  *lon* and *lat* columns.

- algorithm:

  Climatch method algorithm selected from "euclidean" (default) or
  "closest_standard_score".

- d_max:

  Maximum range distance, in kilometres, used when matching occurrence
  points to nearest climate data points/cells (default = 50).

- sd_data:

  Optional `data.frame` for calculating the standard deviation for
  climate variable, or a `vector` of pre-calculated values.

- as_score:

  Logical to indicate whether to generate a score 0-10 (default = TRUE)
  or values 0-1 (FALSE).

- ...:

  Additional parameters.

## Value

A "Climatch" model S4 object containing slots:

- `method`:

  SDM method: "climatch".

- `algorithm`:

  Algorithm: "euclidean" or "closest_standard_score").

- `variables`:

  List of climate (or environmental) variable names.

- `sd`:

  The standard deviation of each variable calculated via the climate
  data (*x*) or the *sd_data* when provided.

- `presence`:

  The selected (nearest within range) climate data for each occurrence
  point.

- `coordinates`:

  The coordinates for the selected climate data.

- `as_score`:

  Indication of whether to generate a score 0-10 or values 0-1.

## References

ABARES (2020). Climatch v2.0 User Manual. Canberra.
<https://climatch.cp1.agriculture.gov.au/> Accessed: November 2021.
