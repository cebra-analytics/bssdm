# Calculate Multivariate Environmental Similarity

Calculate Multivariate Environmental Similarity and most
dissimilar/similar variables with respect to a reference dataset, for a
set of environmental variables.

## Usage

``` r
mess(x, ref, full = FALSE, filename = "", ...)

# S3 method for class 'Raster'
mess(x, ref, full = FALSE, filename = "", ...)

# S3 method for class 'SpatRaster'
mess(x, ref, full = FALSE, filename = "", ...)

# S3 method for class 'data.frame'
mess(x, ref, full = FALSE, ...)

# S3 method for class 'matrix'
mess(x, ref, full = FALSE, ...)
```

## Arguments

- x:

  Climate (or environmental) data as a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  `raster::Raster*`, `list`, `data.frame`, or `matrix` where each
  layer/element/column represents focal values of an environmental
  variable.

- ref:

  A `data.frame`, `matrix`, or `list` where each column/element
  represents reference values for an environmental variable
  (corresponding to those given in `x`).

- full:

  Logical to indicate whether similarity values should be returned for
  all variables. If `FALSE` (the default), then only the minimum
  similarity scores across variables will be returned.

- filename:

  Optional filename for writing spatial raster output (i.e., the `mess`
  raster only). Default is "".

- ...:

  Additional parameters.

## Value

If `x` is a `SpatRaster` or `Raster*` object, this function returns a
list containing:

- `mess`: a `SpatRaster` layer giving the minimum similarity value
  across all variables for each location (i.e. the MESS);

- `mess_by_variable`: a `SpatRaster` giving the environmental
  similarities for each variable in `x` (only included when
  `full=TRUE`);

- `mod`: a factor `SpatRaster` layer indicating which variable was most
  dissimilar to its reference range (i.e. the MoD map, Elith et al.
  2010); and

- `mos`: a factor `SpatRaster` layer indicating which variable was most
  similar to its reference range.

If `x` is a `list`, `matrix` or `data.frame`, the function will return a
list as above, but with multilayer and single layer `SpatRaster` objects
replaced by matrix and vector objects, respectively.

## Details

`mess` uses the MESS algorithm described in Appendix S3 of Elith et al.
2010.

## References

Elith, J., Kearney, M., and Phillips, S. (2010). The art of modelling
range-shifting species. *Methods in Ecology and Evolution*, 1: 330-342.
[doi:10.1111/j.2041-210X.2010.00036.x](https://doi.org/10.1111/j.2041-210X.2010.00036.x)

## Examples

``` r
if (FALSE) { # \dontrun{
library(geodata)
library(terra)
bio <- worldclim_global("bio", res = 10, path = tempdir())
aus <- gadm("AUS", level = 0, resolution = 2, path = tempdir())
occ <- spatSample(aus, size = 100, method = "random")
ref <- terra::extract(bio, occ, ID = FALSE)
m <- mess(bio, ref, full = TRUE)

# Plot outputs
plot(m)
plot(m, which = "mod")
plot(m, which = "mos")
} # }
```
