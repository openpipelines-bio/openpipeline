import mudata
import anndata
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from mudata import MuData
from pathlib import Path
from pandas.testing import assert_frame_equal
from typing import Literal
from .typing import AnnotationObjectOrPathLike


def _read_if_needed(anndata_mudata_path_or_obj):
    if isinstance(anndata_mudata_path_or_obj, (str, Path)):
        return mudata.read(anndata_mudata_path_or_obj)
    if isinstance(anndata_mudata_path_or_obj, (mudata.MuData, anndata.AnnData)):
        return anndata_mudata_path_or_obj.copy()
    raise AssertionError("Expected 'Path', 'str' to MuData/AnnData "
                         "file or MuData/AnnData object.")

def _assert_same_annotation_object_class(left, right):
    assert type(left) == type(right), (f"Two objects are not of the same class:"
                                       f"\n[Left]:{type(left)}\n[right]:{type(right)}")


def assert_mudata_modality_keys_equal(left, right):
    left_keys = set(left.mod.keys())
    right_keys = set(right.mod.keys())
    if left_keys!= right_keys:
        raise AssertionError("MuData modalities differ:"
                             f"\n[left]:{left_keys}\n[right]:{right_keys}")

def assert_shape_equal(left: AnnotationObjectOrPathLike, right: AnnotationObjectOrPathLike):
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    if left.shape != right.shape:
        raise AssertionError(f"{type(left).__name__} shapes differ:"
                             f"\n[left]:{left.shape}\n[right]:{right.shape}")
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            assert_shape_equal(modality, right[mod_name])
 

def assert_obs_names_equal(left: AnnotationObjectOrPathLike, right: AnnotationObjectOrPathLike, 
                           *args, **kwargs):
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    pd.testing.assert_index_equal(left.obs_names, right.obs_names, *args, **kwargs)
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            assert_obs_names_equal(modality, right[mod_name])


def assert_var_names_equal(left: AnnotationObjectOrPathLike, right: AnnotationObjectOrPathLike, 
                           *args, **kwargs):
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    pd.testing.assert_index_equal(left.var_names, right.var_names, *args, **kwargs)
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            assert_var_names_equal(modality, right[mod_name])


def assert_annotation_frame_equal(annotation_attr: Literal["obs", "var"], 
                                   left: AnnotationObjectOrPathLike, right: AnnotationObjectOrPathLike, 
                                   sort=False, *args, **kwargs):
    if not annotation_attr in ("obs", "var"):
        raise ValueError("annotation_attr should be 'obs', or 'var'")
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    left_frame, right_frame = getattr(left, annotation_attr), getattr(right, annotation_attr)
    if sort:
        left_frame, right_frame = left_frame.sort_index(inplace=False), right_frame.sort_index(inplace=False)
    assert_frame_equal(left_frame, right_frame, *args, **kwargs)
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            assert_annotation_frame_equal(annotation_attr, modality, 
                                           right[mod_name], sort=sort, *args, **kwargs)

def _assert_layer_equal(left, right):
    if issparse(left):
        if not issparse(right):
           raise AssertionError("Layers differ:\n[left]: sparse\n[right]: not sparse")
        if left.getformat() != right.getformat():
            raise AssertionError("Layers format differ:"
                                 f"\n[left]:{left.getformat()}\n[right]: {right.getformat()}")
        assert np.all(left.indices == right.indices), "Layers differ: indices are not the same"
        assert np.all(left.indptr == right.indptr), "Layers differ: index pointers are not the same"
        np.testing.assert_allclose(left.data, right.data, 
                                  err_msg="Layers data differs.", equal_nan=True)
    else:
        if issparse(right):
            raise AssertionError("Layers differ:\n[left]: not sparse\n[right]: sparse")
        np.testing.assert_allclose(left, right, 
                                   err_msg="Layers data differs.", equal_nan=True)
        

def assert_layers_equal(left: AnnotationObjectOrPathLike,
                        right: AnnotationObjectOrPathLike):
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    if left.raw is not None:
       _assert_layer_equal(left.raw, right.raw)
    else:
        if right.raw:
            raise AssertionError("Layer .raw differs: "
                                 f"\n[left]:{left.raw}\n[right]:{right}")
    if left.X is not None:
        _assert_layer_equal(left.X, right.X)
    if left.layers:
        assert right.layers and (left.layers.keys() == right.layers.keys()), \
        "Avaiable layers differ:" \
        f"\n[left]:{left.layers}\n[right]{right.layers}"
        for layer_name, layer in left.layers.items():
            _assert_layer_equal(layer, right.layers[layer_name])
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            assert_layers_equal(modality, right[mod_name])


def assert_annotation_objects_equal(left: AnnotationObjectOrPathLike,
                                    right: AnnotationObjectOrPathLike,
                                    check_data=True):
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    assert_shape_equal(left, right)
    assert_annotation_frame_equal("obs", left, right)
    assert_annotation_frame_equal("var", left, right)
    if check_data:
        assert_layers_equal(left, right)
