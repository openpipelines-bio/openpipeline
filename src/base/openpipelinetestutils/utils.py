from .typing import AnnotationObject
from typing import Union, Literal
from operator import attrgetter
from anndata import AnnData
from mudata import MuData
from itertools import product


def remove_annotation_column(annotation_object: AnnotationObject,
                             column_names: list[str] | str,
                             axis: Union[Literal["obs"], Literal["var"], 0, 1],
                             modality_name: str | None = None):
    if isinstance(annotation_object, AnnData) and modality_name is not None:
        raise ValueError("Cannot specify modality when object is of type AnnData.")
    if isinstance(column_names, str):
        column_names = [str(column_names)] # str to make a copy
    axis_strings = {
        "var": "var",
        "obs": "obs",
        0: "obs",
        1: "var"
    }
    axis_string = axis_strings[axis]
    axis_getter = attrgetter(axis_string)
    axis_setter = lambda obj, value: setattr(obj, axis_string, value)
    if not modality_name:
        axis_setter(annotation_object, axis_getter(annotation_object).drop(column_names,
                                                                        axis="columns",
                                                                        inplace=False))
    if isinstance(annotation_object, MuData):
        modality_names = [modality_name] if modality_name else list(annotation_object.mod.keys())
        extra_cols_to_remove = [f"{mod_name}:{column_name}" for mod_name, column_name
                                in product(modality_names, column_names)]
        axis_setter(annotation_object, axis_getter(annotation_object).drop(extra_cols_to_remove,
                                                                           axis="columns",
                                                                           inplace=False))

        for mod_name in modality_names:
            modality = annotation_object.mod[mod_name]
            new_modality = remove_annotation_column(modality, column_names, 
                                                    axis=axis, modality_name=None)
            annotation_object.mod[mod_name] = new_modality 
    return annotation_object
