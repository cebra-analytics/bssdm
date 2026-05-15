# S4 class to represent a climatch SDM.

The model object (S4) class component for an implementation of the
ABARES Climatch species distribution modelling (SDM) method (ABARES,
2020).

## Slots

- `method`:

  SDM method: "climatch".

- `algorithm`:

  Algorithm: "euclidean" or "closest_standard_score".

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
