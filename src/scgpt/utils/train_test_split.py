import anndata as ad
import numpy as np
from sklearn.model_selection import train_test_split
from scipy.sparse import issparse


def extract_data_from_anndata(
    adata: ad.AnnData,
    cell_type_key: str = "celltype",
    batch_id_key: str = "batch_id",
    input_layer_key: str = "X_binned"
) -> dict:
    """
    Extracts data from an AnnData object and returns a dictionary containing the extracted data.

    Parameters:
        adata (ad.AnnData): The AnnData object containing the data.
        cell_type_key (str): The key for the cell type information in the AnnData object. Default is "celltype".
        batch_id_key (str): The key for the batch ID information in the AnnData object. Default is "batch_id".
        input_layer_key (str): The key for the input layer in the AnnData object. Default is "X_binned".

    Returns:
        dict: A dictionary containing the extracted data with the following keys:
            - "all_counts": The input layer data.
            - "celltypes_labels": The cell type labels.
            - "batch_ids": The batch IDs.
    """
    
    all_counts = (
        adata.layers[input_layer_key].A
        if issparse(adata.layers[input_layer_key])
        else adata.layers[input_layer_key]
    )

    celltypes_labels = adata.obs[cell_type_key].tolist()
    celltypes_labels = np.array(celltypes_labels)

    batch_ids = adata.obs[batch_id_key].tolist()
    batch_ids = np.array(batch_ids)
    
    return {
        "all_counts": all_counts,
        "celltypes_labels": celltypes_labels,
        "batch_ids": batch_ids
    }


def split_data_train_test(
    adata: ad.AnnData,
    test_size: float,
    shuffle: bool = True,
) -> dict:
    """
    Split the data into training and testing sets.

    Parameters:
        adata (ad.AnnData): The input AnnData object containing the data.
        test_size (float): The proportion of the data to include in the test set.
        shuffle (bool, optional): Whether to shuffle the data before splitting. Defaults to True.

    Returns:
        dict: A dictionary containing the train and test data, cell type labels, and batch labels.
    """

    data_dict = extract_data_from_anndata(adata)
  
    (
    train_data,
    test_data,
    train_celltype_labels,
    test_celltype_labels,
    train_batch_labels,
    test_batch_labels
    ) = train_test_split(
        data_dict["all_counts"],
        data_dict["celltypes_labels"],
        data_dict["batch_ids"],
        test_size=test_size,
        shuffle=shuffle
    )

    return {
        "train_data": train_data,
        "test_data": test_data,
        "train_celltype_labels": train_celltype_labels,
        "test_celltype_labels": test_celltype_labels,
        "train_batch_labels": train_batch_labels,
        "test_batch_labels": test_batch_labels,
    }
    

def split_data_train_test_val(
    adata: ad.AnnData,
    test_size: float,
    val_size: float,
    shuffle: bool = True,
) -> dict:
    """
    Split the data into training, testing, and validation sets.

    Parameters:
        adata (ad.AnnData): The input AnnData object containing the data.
        test_size (float): The proportion of the data to include in the test set.
        val_size (float): The proportion of the data to include in the validation set.
        shuffle (bool, optional): Whether to shuffle the data before splitting. Defaults to True.

    Returns:
        dict: A dictionary containing the split data and labels.

    """
    data_dict = extract_data_from_anndata(adata)
  
    (
    train_data,
    test_data,
    train_celltype_labels,
    test_celltype_labels,
    train_batch_labels,
    test_batch_labels
    ) = train_test_split(
        data_dict["all_counts"],
        data_dict["celltypes_labels"],
        data_dict["batch_ids"],
        test_size=test_size,
        shuffle=shuffle
    )

    (
    train_data,
    val_data,
    train_celltype_labels,
    val_celltype_labels,
    train_batch_labels,
    val_batch_labels
    ) = train_test_split(
        train_data,
        train_celltype_labels,
        train_batch_labels,
        test_size=val_size/(1-test_size),
        shuffle=shuffle
    )

    return {
        "train_data": train_data,
        "test_data": test_data,
        "val_data": val_data,
        "train_celltype_labels": train_celltype_labels,
        "test_celltype_labels": test_celltype_labels,
        "val_celltype_labels": val_celltype_labels,
        "train_batch_labels": train_batch_labels,
        "test_batch_labels": test_batch_labels,
        "val_batch_labels": val_batch_labels,
    }


