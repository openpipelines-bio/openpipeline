# KNN Labels Transfer

This component transfers labels from a reference dataset to an query using a weighted K-Nearest Neighbors classification approach based on aproximate neighbors search of [PyNNDescent](https://pynndescent.readthedocs.io/en/latest/how_to_use_pynndescent.html).

# Input

- `.h5mu` files for reference and query as specified in [config](./config.vsh.yaml)
- List of targets (labels) to predict from `.obs` slot of the reference. Common choice for HLCA would be `["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"]`
- Number of nearest neighbors to use for prediction. We recommend using 30-50 neighbors for HLCA size dataset

For the description of other input arguments, refer to the [config](./config.vsh.yaml).

# Output

`.h5mu` file with the query dataset that contains the transferred annotations
- The annotations for each of the `targets` will be stored as new columns in the `.obs` data frame with the suffix `_pred` by default
- The uncertainty of the prediction for each target will be stored as another column in `.obs` with the suffix `_uncertainty`
- `labels_transfer` entry is added to the `.uns` dictionary of the output file, which contains information about the method used and its parameters. Keys of this dictionary are labels from `targets` list. Values are dictionaries with the following keys:

| Key         | Description                   |
|-------------|------------------------|
| method      | Always the string "KNN_pynndescent" |
| n_neighbors | Integer indicating the number of neighbors used for the prediction |
| reference   | File or a link used as `reference` input argument |