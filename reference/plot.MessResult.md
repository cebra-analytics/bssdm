# Plot method for MessResult objects

Produces a tmap plot for outputs of
[`mess`](https://cebra-analytics.github.io/bssdm/reference/mess.md).

## Usage

``` r
# S3 method for class 'MessResult'
plot(x, which = "mess", ...)
```

## Arguments

- x:

  A `MessResult` object returned by
  [`mess`](https://cebra-analytics.github.io/bssdm/reference/mess.md).

- which:

  Character string specifying which output to plot. One of `"mess"`
  (default), `"mod"` (most dissimilar variable), or `"mos"` (most
  similar variable).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`. Called for its side-effect of rendering a tmap
plot.

## Note

Requires the tmap package to be installed.
