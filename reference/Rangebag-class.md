# S4 class to represent a range bagging SDM.

The model object (S4) class component for an implementation of the range
bagging species distribution modelling (SDM) method (Drake, 2015).

## Slots

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
