# XGBoost Labels Transfer

This component transfers labels from a reference dataset to an query using a boosting method [XGBoost](https://xgboost.readthedocs.io/en/stable/).

# Input

- `.h5mu` files for reference and query as specified in [config](./config.vsh.yaml)
- List of targets (labels) to predict from `.obs` slot of the reference. Common choice for HLCA would be `["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"]`
- models output directory

For the description of other input arguments, refer to the [config](./config.vsh.yaml).

# Output

`.h5mu` file with the query dataset that contains the transferred annotations
- The annotations for each of the `targets` will be stored as new columns in the `.obs` data frame with the suffix `_pred` by default
- The uncertainty of the prediction for each target will be stored as another column in `.obs` with the suffix `_uncertainty`
- `labels_transfer` entry is added to the `.uns` dictionary of the output file, which contains information about the method used and its parameters. Keys of this dictionary are labels from `targets` list. Values are dictionaries with the following keys:

| Key         | Description                   |
|-------------|------------------------|
| method      | Always the string "XGBClassifier" |
| learning_rate | Training parameter of XGBoost classifier used |
| min_split_loss   | Training parameter of XGBoost classifier used |
| max_depth   | Training parameter of XGBoost classifier used |
| min_child_weight   | Training parameter of XGBoost classifier used |
| max_delta_step   | Training parameter of XGBoost classifier used|
| subsample   | Training parameter of XGBoost classifier used |
| sampling_method   | Training parameter of XGBoost classifier used |
| colsample_bytree   |  Training parameter of XGBoost classifier used |
| colsample_bylevel   | Training parameter of XGBoost classifier used |
| colsample_bynode   | Training parameter of XGBoost classifier used|
| reg_lambda   | Training parameter of XGBoost classifier used |
| reg_alpha   | Training parameter of XGBoost classifier used |
| scale_pos_weight   | Training parameter of XGBoost classifier used |

- Models saved in models output directory (`./model` by default). File for each target is named `classifier_<target>.xgb`
- `model_info.json` file saved in models output directory. It has the following format:
1. The first level key is always the string "classifier_info". This is done for the possibility of further development
2. The second level keys are labels from the list of `targets`. The values are dictionaries with the following format:

| Key | Description | Example |
|-------------|------------------------|----------|
| filename      | Name of the file containing classifier in models output directory| classifier_ann_level_1.xgb |
| labels      | List of unique labels in the corresponding column in reference | ["Stroma", "Endothelial", "Immune"] |
| obs_column      | Output column for the predicted labels in the output dataset | ann_level_1_pred |
| model_params      | Dictionary with the training parameters, the same that in `.uns` dictionary of the output dataset | {"learning_rate": 0.3, ...} | 

# Notes

- If there are suitable files with classifiers in the models output directory, the trained models from these files would be used for the labels transfer
- New models are trained only for targets, for which trained classifiers are not found
- If you want to retrain all the classifier, use the `force_retrain` parameter during running of the component. Be carefut: it may overwrite the files in models output directory
