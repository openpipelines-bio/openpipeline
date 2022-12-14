import json
import os
from typing import Optional

import mudata
import numpy as np
import scanpy as sc
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.preprocessing import LabelEncoder


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "reference": "https://zenodo.org/record/6337966/files/HLCA_emb_and_metadata.h5ad",
    "targets": ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"],
    "modality": "rna",
    "reference_obsm_key": "X_integrated_scanvi",
    "query_obsm_key": "X_integrated_scanvi",
    "output": "foo.h5mu",
    "model_output": "model",
    "obs_output_suffix": "_pred",
    "force_retrain": False,
    "use_gpu": True,
    "verbosity": 1,
    "learning_rate": 0.3,
    "min_split_loss": 0,
    "max_depth": 6,
    "min_child_weight": 1,
    "max_delta_step": 0,
    "subsample": 1,
    "sampling_method": "uniform",
    "colsample_bytree": 1,
    "colsample_bylevel": 1,
    "colsample_bynode": 1,
    "reg_lambda": 1,
    "reg_alpha": 0,
    "scale_pos_weight": 1,
}
### VIASH END

training_params = {
    "learning_rate": par["learning_rate"],
    "min_split_loss": par["min_split_loss"],
    "max_depth": par["max_depth"],
    "min_child_weight": par["min_child_weight"],
    "max_delta_step": par["max_delta_step"],
    "subsample": par["subsample"],
    "sampling_method": par["sampling_method"],
    "colsample_bytree": par["colsample_bytree"],
    "colsample_bylevel": par["colsample_bylevel"],
    "colsample_bynode": par["colsample_bynode"],
    "reg_lambda": par["reg_lambda"],
    "reg_alpha": par["reg_alpha"],
    "scale_pos_weight": par["scale_pos_weight"],
}

def encode_labels(y):
    labels_encoder = LabelEncoder()
    labels_encoder.fit(y)
    
    return labels_encoder.transform(y), labels_encoder


def get_model_eval(xgb_model, X_test, y_test, labels_encoder):
    
    preds = xgb_model.predict(X_test)
    
    cr = classification_report(labels_encoder.inverse_transform(y_test),
                               labels_encoder.inverse_transform(preds),
                               output_dict=True)
    cr_df = pd.DataFrame(cr).transpose()

    return cr_df


def train_test_split_adata(adata, labels):
    
    train_data = pd.DataFrame(data=adata.X, index=adata.obs_names)

    X_train, X_test, y_train, y_test = train_test_split(
        train_data, labels, test_size=0.2, random_state=42, stratify=labels)
    
    return X_train, X_test, y_train, y_test


def train_xgb_model(X_train, y_train, gpu=True) -> xgb.XGBClassifier:
    
    n_classes = len(np.unique(y_train))
    objective = "binary:logistic" if n_classes == 2 else "multi:softprob"
    
    tree_method = "gpu_hist" if gpu else "hist"        
    xgbc = xgb.XGBClassifier(tree_method=tree_method, objective=objective, **training_params)
    xgbc.fit(X_train, y_train)
    
    return xgbc


def build_classifier(X, y, labels_encoder, label_key, eval_verbosity: Optional[int] = 1, gpu=True) -> xgb.XGBClassifier:
    # Adata prep
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    #Note: Do we need a new train-test split for each classifier?
    
    # Model training
    xgb_model = train_xgb_model(X_train, y_train, gpu=gpu)
    
    # Model eval
    if eval_verbosity != 0:
        cr_df = get_model_eval(xgb_model, X_test, y_test, labels_encoder)
        
        if eval_verbosity == 2:
            print(cr_df)

        else:
            overall_accuracy = cr_df["support"]["accuracy"]
            low_prec_key = cr_df.precision.idxmin()
            low_prec_val = cr_df.precision.min()
            low_rec_key = cr_df.recall.idxmin()
            low_rec_val = cr_df.recall.min()
            low_f1_key = cr_df["f1-score"].idxmin()
            low_f1_val = cr_df["f1-score"].min()

            print("")
            print(f"Summary stats for {label_key} model:")
            print(f"Overall accuracy: {overall_accuracy}")
            print(f"Min. precision: {low_prec_key}: {low_prec_val}")
            print(f"Min. Recall: {low_rec_key}: {low_rec_val}")
            print(f"Min. F1-score: {low_f1_key}: {low_f1_val}")
            print("")
            
    return xgb_model


def build_ref_classifiers(adata_reference, targets, model_path,
                          eval_verbosity: Optional[int] = 1, gpu: Optional[bool] = True) -> None:
    """
    This function builds xgboost classifiers on a reference embedding for a designated number of 
    adata_reference.obs columns. Classifier .xgb files and a model_info.json file is written to the `model_path`
    directory. Model evaluation is printed to stdout.
    
    Inputs:
        * `adata_reference`: The AnnData object that was used to train the reference model
        * `model_path`: The reference model directory where the classifiers will also be stored
        * `eval_verbosity`: The verbosity level for evaluation of the classifier from the range [0;2].
        * `gpu`: Boolean indicating whether a gpu is available for classifier training
        
        
    Example:
    ```
    >>> adata
    AnnData object with n_obs × n_vars = 700 × 765
    obs: "ann_finest_level", "ann_level_1"
    
    >>> os.listdir("/path/to/model")
    model_params.pt*
    
    >>> build_ref_classifiers(adata, "path/to/model", eval_verbosity=1, gpu=True)
    >>> os.listdir("/path/to/model")
    classifier_ann_finest_level.xgb*    model_info.json*
    classifier_ann_level_1.xgb*         model_params.pt* 
    ```
    
    """
    # Check inputs
    if not isinstance(eval_verbosity, int):
        raise TypeError("`eval_verbosity` should be an integer between 0 and 2.")
        
    if eval_verbosity < 0 or eval_verbosity > 2:
        raise ValueError("`eval_verbosity` should be an integer between 0 and 2.")

    if par["reference_obsm_key"] is None:
        X_reference = adata_reference.X
    else:
        X_reference = adata_reference.obsm[par["reference_obsm_key"]]
    
    # Map from name of classifier to file names
    classifiers = dict()
    
    for label in targets:
        if label not in adata_reference.obs:
            raise ValueError(f"{label} is not in the `adata` object passed!")

        filename = "classifier_" + label + ".xgb"
        
        labels, labels_encoder = encode_labels(adata_reference.obs[label])
        print("Classes:", labels_encoder.classes_)
        
        print(f"Building classifier for {label}...")
        xgb_model = build_classifier(
            X=X_reference,
            y=labels,
            labels_encoder=labels_encoder,
            label_key=label,
            eval_verbosity=eval_verbosity,
            gpu=gpu
        )
        
        # Save classifier
        xgb_model.save_model(os.path.join(model_path, filename))
        
        # Store classifier info
        classifiers[label] = {
            "filename": filename,
            "labels": labels_encoder.classes_.tolist()
        }
        
    # Store model_info.json file
    model_info = {
        "model_params": training_params,
        "classifier_info": classifiers
    }
    
    # Write file
    json_string = json.dumps(model_info, indent=4)
    
    with open(model_path + "/model_info.json", "w") as f:
        f.write(json_string)  # TODO: if there is already the existing file, it will overwrite the whole info


def main():
    mdata = mudata.read(par["input"].strip())
    adata = mdata.mod[par["modality"]]

    adata_reference = sc.read(par["reference"], backup_url=par["reference"])

    # If classifiers for targets are in the model_output directory, simply open them and run (unless `retrain` != True)
    # If some classifiers are missing, train and save them first
    # Predict and save the query data

    if par["query_obsm_key"] is None:
        query = adata.X
    else:
        query = adata.obsm[par["query_obsm_key"]]

    targets_to_train = []

    for target in par["targets"]:
        if f"classifier_{target}.xgb" not in os.listdir(par["model_output"]) or par["force_retrain"]:
            targets_to_train.append(target)

    print("Verbosity:", par["verbosity"], type(par["verbosity"]))
    # TODO: align features names?
    build_ref_classifiers(adata_reference, targets_to_train, model_path=par["model_output"], gpu=par["use_gpu"], eval_verbosity=par["verbosity"])

    if "labels_transfer" not in adata.uns:
        adata.uns["labels_transfer"] = {}

    for target in par["targets"]:
        predicted_label_col_name = target + par["obs_output_suffix"]



        # adata.obs[predicted_label_col_name] = prediction
        # adata.obs[target + "_uncertainty"] = uncertainty
        
        adata.uns["labels_transfer"][predicted_label_col_name] = {
            "method": "XGBClassifier",
            **training_params  # TODO: if classifier was trained previously, this might corrupt its info
        }

    mdata.mod[par['modality']] = adata
    mdata.update()
    mdata.write_h5mu(par['output'].strip())

if __name__ == "__main__":
    main()
