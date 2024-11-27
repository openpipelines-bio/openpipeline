import mudata
import anndata
import pandas as pd
import numpy as np
from scipy.sparse import issparse, spmatrix
from mudata import MuData
from pathlib import Path
from pandas.testing import assert_frame_equal
from typing import Literal
from .typing import AnnotationObjectOrPathLike
from functools import singledispatch


def _read_if_needed(anndata_mudata_path_or_obj):
    if isinstance(anndata_mudata_path_or_obj, (str, Path)):
        return mudata.read(str(anndata_mudata_path_or_obj)) # TODO: remove when mudata fixes PAth bug
    if isinstance(anndata_mudata_path_or_obj, (mudata.MuData, anndata.AnnData)):
        return anndata_mudata_path_or_obj.copy()
    raise AssertionError("Expected 'Path', 'str' to MuData/AnnData "
                         "file or MuData/AnnData object.")

def _assert_same_annotation_object_class(left, right):
    assert type(left) == type(right), (f"Two objects are not of the same class:"
                                       f"\n[Left]:{type(left)}\n[right]:{type(right)}")
    
def _promote_dtypes(left, right):
    # Create new DataFrames to avoid modifying the original ones
    left_aligned = left.copy()
    right_aligned = right.copy()
    
    for column in left.columns:
        l_dtype = left[column].dtype
        r_dtype = right[column].dtype
        
        if l_dtype == r_dtype:
            # No need to modify dtypes that are already the same
            continue
        if not all(map(pd.api.types.is_any_real_numeric_dtype, (r_dtype, l_dtype))):
            # Do not try casting without dtypes that do not represent real numbers
            continue
        is_extension = pd.api.types.is_extension_array_dtype(l_dtype)
        if is_extension and not pd.api.types.is_extension_array_dtype(r_dtype):
            continue
        numpy_dtype_l = l_dtype.type if is_extension else l_dtype
        numpy_dtype_r = r_dtype.type if is_extension else r_dtype 
        # At this point we should have only integer or float dtypes 
        common_dtype = np.promote_types(numpy_dtype_l, numpy_dtype_r)
        if is_extension:
            left_aligned[column] = pd.array(left[column], dtype=common_dtype)
            right_aligned[column] = pd.array(right[column], dtype=common_dtype)
        else:
            left_aligned[column] = left[column].astype(common_dtype)
            right_aligned[column] = right[column].astype(common_dtype)
    
    return left_aligned, right_aligned


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


def _assert_frame_equal(left, right, sort=False, promote_precicion=False, *args, **kwargs):
    if sort:
        left, right = left.sort_index(inplace=False), right.sort_index(inplace=False)
        left, right = left.sort_index(axis=1, inplace=False), right.sort_index(axis=1, inplace=False)
        
    if promote_precicion:
        left, right = _promote_dtypes(left, right)
        assert_frame_equal(left, right, check_exact=False, atol=1e-3, *args, **kwargs)
    else:
        assert_frame_equal(left, right, *args, **kwargs)

def assert_annotation_frame_equal(annotation_attr: Literal["obs", "var"], 
                                   left: AnnotationObjectOrPathLike, right: AnnotationObjectOrPathLike, 
                                   sort=False,
                                   promote_precicion=False,
                                   *args, **kwargs):
    if not annotation_attr in ("obs", "var"):
        raise ValueError("annotation_attr should be 'obs', or 'var'")
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    left_frame, right_frame = getattr(left, annotation_attr), getattr(right, annotation_attr)
    _assert_frame_equal(left_frame, right_frame, sort=sort, promote_precicion=promote_precicion, *args, **kwargs)
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            assert_annotation_frame_equal(annotation_attr, modality, 
                                           right[mod_name], sort=sort, promote_precicion=promote_precicion, *args, **kwargs)

def _assert_layer_equal(left, right):
    if issparse(left):
        if not issparse(right):
           raise AssertionError("Layers differ:\n[left]: sparse\n[right]: not sparse")
        if left.getformat() != right.getformat():
            raise AssertionError("Layers format differ:"
                                 f"\n[left]:{left.getformat()}\n[right]: {right.getformat()}")
        assert np.all(left.indices == right.indices), "Layers differ: indices are not the same"
        assert (left.indices.dtype == left.indptr.dtype), \
            "Problem with layer from 'left': dtype of sparce array 'inptr' and 'indices' are not the same"
        assert (right.indices.dtype == right.indptr.dtype), \
            "Problem with layer from 'right': dtype of sparce array 'inptr' and 'indices' are not the same"
        assert (left.indices.dtype == right.indices.dtype), \
            "Layers differ: dtype of sparce matrix 'indices' are different."
        assert (left.indptr.dtype == right.indptr.dtype), \
            "Layers differ: dtype of sparce matrix 'indptr' are different."
        assert np.all(left.indptr == right.indptr), "Layers differ: index pointers are not the same"
        np.testing.assert_allclose(left.data, right.data, rtol=1e-5,
                                  err_msg="Layers data differs.", equal_nan=True)
    else:
        if issparse(right):
            raise AssertionError("Layers differ:\n[left]: not sparse\n[right]: sparse")
        np.testing.assert_allclose(left, right, 
                                   rtol=1e-5,
                                   err_msg="Layers data differs.",
                                   equal_nan=True)
        

def assert_layers_equal(left: AnnotationObjectOrPathLike,
                        right: AnnotationObjectOrPathLike):
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    if left.raw is not None:
        try:
            _assert_layer_equal(left.raw, right.raw)
        except AssertionError as e:
            e.add_note(".raw is different")
            raise
    else:
        if right.raw:
            raise AssertionError("Layer .raw differs: "
                                 f"\n[left]:{left.raw}\n[right]:{right}")
    if left.X is not None:
        try:
            _assert_layer_equal(left.X, right.X)
        except AssertionError as e:
            e.add_note("X is different.")
            raise
    if left.layers:
        assert right.layers and (left.layers.keys() == right.layers.keys()), \
        "Avaiable layers differ:" \
        f"\n[left]:{left.layers}\n[right]{right.layers}"
        for layer_name, layer in left.layers.items():
            try:
                _assert_layer_equal(layer, right.layers[layer_name])
            except AssertionError as e:
                e.add_note(f"Layer {layer_name} is different")
                raise
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            assert_layers_equal(modality, right[mod_name])



def assert_multidimensional_annotation_equal(annotation_attr: Literal["obsm", "varm"],
                                             left, right, sort=False):
    if not annotation_attr in ("obsm", "varm"):
        raise ValueError("annotation_attr should be 'obsm', or 'varm'")
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)

    @singledispatch
    def _assert_multidimensional_value_equal(left, right, **kwargs):
        raise NotImplementedError("Unregistered type found while asserting")
    
    @_assert_multidimensional_value_equal.register
    def _(left: pd.DataFrame, right, **kwargs):
        _assert_frame_equal(left, right, **kwargs)
   
    @_assert_multidimensional_value_equal.register(np.ndarray)
    @_assert_multidimensional_value_equal.register(spmatrix)
    def _(left, right, **kwargs):
        # Cannot sort sparse and dense matrices so ignore sort param
        _assert_layer_equal(left, right)

    left_dict, right_dict = getattr(left, annotation_attr), getattr(right, annotation_attr)
    left_keys, right_keys = left_dict.keys(), right_dict.keys()
    assert left_keys == right_keys, f"Keys of {annotation_attr} differ:\n[left]:{left_keys}\n[right]:{right_keys}"
    for left_key, left_value in left_dict.items():
        try:
            _assert_multidimensional_value_equal(left_value, right_dict[left_key], sort=sort)
        except AssertionError as e:
            e.add_note(f"Failing key: {left_key}")
            raise
    if isinstance(left, MuData):
        assert_mudata_modality_keys_equal(left, right)
        for mod_name, modality in left.mod.items(): 
            try:
                assert_multidimensional_annotation_equal(annotation_attr ,modality, right[mod_name], sort=sort)
            except AssertionError as e:
                e.add_note(f"Failing modality: {mod_name}")
                raise

def assert_annotation_objects_equal(left: AnnotationObjectOrPathLike,
                                    right: AnnotationObjectOrPathLike,
                                    check_data=True,
                                    sort=True,
                                    promote_precision=False):
    left, right = _read_if_needed(left), _read_if_needed(right)
    _assert_same_annotation_object_class(left, right)
    assert_shape_equal(left, right)
    assert_annotation_frame_equal("obs", left, right, sort=sort, promote_precicion=promote_precision)
    assert_annotation_frame_equal("var", left, right, sort=sort, promote_precicion=promote_precision)
    for slot in ("varm", "obsm"):
        try:
            assert_multidimensional_annotation_equal(slot, left, right, sort=sort)
        except AssertionError as e:
            e.add_note(f"Failing multidimensional slot: {slot}")
            raise
    if check_data:
        assert_layers_equal(left, right)