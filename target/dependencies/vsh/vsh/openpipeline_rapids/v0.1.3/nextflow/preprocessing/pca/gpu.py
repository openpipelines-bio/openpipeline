from contextlib import contextmanager
import numpy as np
import rapids_singlecell as rsc


def _ensure_gpu_compatible_dtype(matrix):
    """Cast integer matrices to float32 so they can be moved to the GPU.

    cupy sparse matrices only support bool/float/complex dtypes, so integer
    count matrices (e.g. raw counts) cannot be transferred as-is. The
    rapids-singlecell operations that run on the transferred data produce float
    values anyway, so casting integer counts to float loses no information.
    """
    if matrix is not None and np.issubdtype(matrix.dtype, np.integer):
        return matrix.astype("float32")
    return matrix


def _cast_transfer_targets(adata, kwargs):
    """Cast the matrices ``anndata_to_GPU`` will transfer to a GPU-safe dtype.

    Mirrors what ``anndata_to_GPU`` moves for the given keyword arguments so
    that matrices it does not touch keep their original dtype.
    """
    if kwargs.get("convert_all"):
        adata.X = _ensure_gpu_compatible_dtype(adata.X)
        for key in list(adata.layers.keys()):
            adata.layers[key] = _ensure_gpu_compatible_dtype(adata.layers[key])
    elif kwargs.get("layer"):
        layer = kwargs["layer"]
        adata.layers[layer] = _ensure_gpu_compatible_dtype(adata.layers[layer])
    else:
        adata.X = _ensure_gpu_compatible_dtype(adata.X)


@contextmanager
def on_gpu(adata, logger=None, *, slots=None, **kwargs):
    """Move an AnnData to the GPU for the duration of the block, then back.

    Extra keyword arguments are forwarded to ``rsc.get.anndata_to_GPU``
    (e.g. ``layer=...`` or ``convert_all=True``).

    ``anndata_to_GPU`` only transfers ``.X`` and ``.layers``; it never touches
    aligned mappings such as ``.obsm`` or ``.varm``. Pass ``slots`` to move
    individual entries of those via ``X_to_GPU``, restoring them with
    ``X_to_CPU`` on exit. ``slots`` maps an AnnData attribute name to the
    key(s) within it, e.g. ``{"obsm": "X_pca"}`` or
    ``{"obsm": ["X_pca", "X_umap"], "varm": ["loadings"]}``. When only ``slots``
    is requested, ``.X``/``.layers`` are left on the CPU.

    Integer matrices among the transferred ``.X``/``.layers`` are cast to
    float32 first, since cupy sparse cannot represent integer dtypes.
    """
    # Components treat an empty --layer as "use .X", and pass it through
    # explicitly rather than omitting the argument. An empty string is not a
    # layer name, so drop it instead of forwarding it to anndata_to_GPU.
    if "layer" in kwargs and not kwargs["layer"]:
        del kwargs["layer"]

    slot_keys = {
        attr: [keys] if isinstance(keys, str) else list(keys)
        for attr, keys in (slots or {}).items()
    }
    transfer_adata = bool(kwargs) or not slot_keys

    def move_slots(transfer_func):
        for attr, keys in slot_keys.items():
            mapping = getattr(adata, attr)
            for key in keys:
                mapping[key] = transfer_func(mapping[key])

    description = _describe_transfer(transfer_adata, kwargs, slot_keys)

    if logger is not None:
        logger.info("Transferring %s to GPU.", description)
    if transfer_adata:
        _cast_transfer_targets(adata, kwargs)
        rsc.get.anndata_to_GPU(adata, **kwargs)
    move_slots(rsc.get.X_to_GPU)

    yield adata

    if logger is not None:
        logger.info("Transferring %s back to CPU.", description)
    if transfer_adata:
        rsc.get.anndata_to_CPU(adata)
    move_slots(rsc.get.X_to_CPU)


def _describe_transfer(transfer_adata, kwargs, slot_keys):
    """Human-readable summary of what ``on_gpu`` moves, for logging."""
    parts = []
    if transfer_adata:
        if kwargs.get("convert_all"):
            parts.append(".X and all layers")
        elif "layer" in kwargs:
            parts.append(f"layer {kwargs['layer']!r}")
        else:
            parts.append(".X and .layers")
    for attr, keys in slot_keys.items():
        keys_str = ", ".join(repr(key) for key in keys)
        parts.append(f".{attr}[{keys_str}]")
    return ", ".join(parts) if parts else "nothing"
