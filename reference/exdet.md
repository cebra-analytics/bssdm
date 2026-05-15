# Extrapolation Detection for Species Distribution Models

Detect and quantify both univariate (Type 1) and multivariate (Type 2)
environmental novelty when projecting species distribution models.

## Usage

``` r
exdet(x, ref, mic = FALSE, filename = "", tol = .Machine$double.eps, ...)

# S3 method for class 'Raster'
exdet(x, ref, mic = FALSE, filename = "", tol = .Machine$double.eps, ...)

# S3 method for class 'SpatRaster'
exdet(x, ref, mic = FALSE, filename = "", tol = .Machine$double.eps, ...)

# S3 method for class 'data.frame'
exdet(x, ref, mic = FALSE, tol = .Machine$double.eps, ...)

# S3 method for class 'matrix'
exdet(x, ref, mic = FALSE, tol = .Machine$double.eps, ...)
```

## Arguments

- x:

  Climate (or environmental) data as a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  `raster::Raster*`, `data.frame`, or `matrix` where each layer/column
  represents focal values of an environmental variable.

- ref:

  A `data.frame`, `matrix`, or `list` where each column/element
  represents reference values for an environmental variable
  (corresponding to those given in `x`).

- mic:

  Logical to indicate whether most influential covariates should be
  returned. If `TRUE`, the function returns which variables are most
  responsible for Type 1 and Type 2 novelty. Default is `FALSE`.

- filename:

  Optional filename for writing spatial raster output (i.e., the `exdet`
  raster only). Default is "".

- tol:

  Tolerance value passed to
  [`mahalanobis()`](https://rdrr.io/r/stats/mahalanobis.html) for matrix
  inversion. Default is `.Machine$double.eps`. See
  [`?solve`](https://rdrr.io/r/base/solve.html) for details.

- ...:

  Additional parameters.

## Value

If `x` is a `SpatRaster` or `Raster*` object, this function returns a
list containing:

- `exdet`: a `SpatRaster` layer giving the extrapolation detection
  scores where values \< 0 indicate univariate novelty (Type 1), values
  between 0 and 1 indicate analog conditions, and values \> 1 indicate
  novel covariate combinations (Type 2);

- `mic1`: a factor `SpatRaster` layer indicating which variable is most
  influential for Type 1 novelty. Value is "Not novel" where no Type 1
  novelty occurs (only included when `mic=TRUE`); and

- `mic2`: a factor `SpatRaster` layer indicating which variable is most
  influential for Type 2 novelty. Value is "Not novel" where no Type 2
  novelty occurs (only included when `mic=TRUE`).

If `x` is a `data.frame` or `matrix`, the function will return a list as
above, but with single layer `SpatRaster` objects replaced by vectors,
and factor rasters replaced by factor vectors.

## Details

`exdet` implements the ExDet (Extrapolation Detection) method described
in Mesgaran et al. (2014). It detects both novel univariate ranges (Type
1 novelty) and novel combinations of covariates (Type 2 novelty) in the
projection area compared to the reference area.

Type 1 novelty is assessed by comparing each environmental variable to
its range in the reference data. Type 2 novelty is assessed using
Mahalanobis distance to detect novel combinations of variables, even
when individual variables are within their reference ranges.

The ExDet score prioritises univariate novelty. If one or more variables
are outside their reference range, Type 2 novelty is not considered.

When `mic=TRUE`, the most influential covariates (MIC) are identified
for both types of novelty. MIC1 identifies which variable contributes
most to Type 1 novelty (most outside its reference range). MIC2
identifies which variable contributes most to Type 2 novelty (via
leave-one-out influence analysis). Areas with no novelty are labeled
"Not novel" in the MIC outputs.

## References

Mesgaran, M. B., Cousens, R. D., and Webber, B. L. (2014). Here be
dragons: a tool for quantifying novelty due to covariate range and
correlation change when projecting species distribution models.
*Diversity and Distributions*, 20(10): 1147-1159.
[doi:10.1111/ddi.12209](https://doi.org/10.1111/ddi.12209)

## Examples

``` r
if (FALSE) { # \dontrun{
library(geodata)
library(terra)
bio <- worldclim_global("bio", res = 10, path = tempdir())
aus <- gadm("AUS", level = 0, resolution = 2, path = tempdir())
occ <- spatSample(aus, size = 100, method = "random")
ref <- terra::extract(bio, occ, ID = FALSE)
ex <- exdet(bio, ref, mic = TRUE)

# Plot outputs
plot(ex)
plot(ex, which = "mic1")
plot(ex, which = "mic2")
} # }
```
