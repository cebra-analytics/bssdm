# Range bagging SDM model building method

The model building component of an implementation of the range bagging
species distribution modelling (SDM) method (Drake, 2015).

## Usage

``` r
rangebag(
  x,
  p,
  n_models = 100,
  n_dim = 2,
  sample_prop = 0.5,
  limit_occur = TRUE,
  ...
)

# S3 method for class 'Raster'
rangebag(
  x,
  p,
  n_models = 100,
  n_dim = 2,
  sample_prop = 0.5,
  limit_occur = TRUE,
  ...
)

# S3 method for class 'SpatRaster'
rangebag(
  x,
  p,
  n_models = 100,
  n_dim = 2,
  sample_prop = 0.5,
  limit_occur = TRUE,
  ...
)

# S3 method for class 'data.frame'
rangebag(
  x,
  p,
  n_models = 100,
  n_dim = 2,
  sample_prop = 0.5,
  limit_occur = TRUE,
  ...
)
```

## Arguments

- x:

  Climate (or environmental) data as a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or `raster::Raster*` (with any CRS), or a `data.frame` with WGS84
  *lon* and *lat* columns.

- p:

  Species occurrence data as a `data.frame` (or `matrix`) with WGS84
  *lon* and *lat* columns. Any points outside the extent of `x` will be
  ignored.

- n_models:

  Number of convex hull models to build in sampled environment space
  (default = 100).

- n_dim:

  Number of dimensions (variables) of sampled convex hull models
  (default = 2).

- sample_prop:

  Proportion of environment data rows sampled for fitting each convex
  hull model (default = 0.5).

- limit_occur:

  Logical to indicate whether to limit occurrence data to one per
  environment data cell (default = TRUE).

- ...:

  Additional parameters.

## Value

A "Rangebag" model S4 object containing slots:

- `method`:

  SDM method: "rangebag".

- `variables`:

  List of climate (or environmental) variable names.

- `presence`:

  The selected climate data corresponding to occurrences points.

- `coordinates`:

  The coordinates for the selected climate data.

- `ch_models`:

  A list of convex hull models (vertices).

## References

Drake, J. M. (2015). Range bagging: a new method for ecological niche
modelling from presence-only data. *Journal of the Royal Society
Interface*, 12(107), 20150086.
[doi:10.1098/rsif.2015.0086](https://doi.org/10.1098/rsif.2015.0086)
