from typing import Union
from mudata import MuData
from anndata import AnnData
from pathlib import Path

AnnotationObject = Union[MuData, AnnData]
AnnotationObjectOrPathLike = Union[AnnotationObject, str, Path]