# Plot method for ExdetResult objects

Produces a tmap plot for outputs of
[`exdet`](https://cebra-analytics.github.io/bssdm/reference/exdet.md).

## Usage

``` r
# S3 method for class 'ExdetResult'
plot(x, which = "exdet", ...)
```

## Arguments

- x:

  An `ExdetResult` object returned by
  [`exdet`](https://cebra-analytics.github.io/bssdm/reference/exdet.md).

- which:

  Character string specifying which output to plot. One of `"exdet"`
  (default), `"mic1"` (most influential covariate for Type 1 novelty),
  or `"mic2"` (most influential covariate for Type 2 novelty).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`. Called for its side-effect of rendering a tmap
plot.

## Note

Requires the tmap package to be installed.
