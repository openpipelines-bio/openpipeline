import sys
from mudata import read_h5ad

## VIASH START
par = {
    "input_source": "source.h5mu",
    "source_modality": "rna",
    "input_target": "target.h5mu",
    "target_modality": None,
    "obs": None,
    "var": None,
    "obsm": None,
    "varm": None,
    "obsp": None,
    "varp": None,
    "uns": None,
    "allow_overwrite": False,
    "output": "output.h5mu",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

target_modality = par["target_modality"] or par["source_modality"]

logger.info(
    "Reading modality '%s' from source file '%s'",
    par["source_modality"],
    par["input_source"],
)
try:
    source_mod = read_h5ad(par["input_source"], mod=par["source_modality"])
except KeyError:
    raise ValueError(
        f"Modality '{par['source_modality']}' does not exist in source file "
        f"'{par['input_source']}'."
    )

logger.info(
    "Reading modality '%s' from target file '%s'",
    target_modality,
    par["input_target"],
)
try:
    target_mod = read_h5ad(par["input_target"], mod=target_modality)
except KeyError:
    raise ValueError(
        f"Modality '{target_modality}' does not exist in target file "
        f"'{par['input_target']}'."
    )

# .obs/.var are DataFrames (column access), .obsm/.varm/.obsp/.varp are array
# containers, and .uns is a dict -- all support key-based get/set via getattr.
_slots = [
    ("obs", par["obs"]),
    ("var", par["var"]),
    ("obsm", par["obsm"]),
    ("varm", par["varm"]),
    ("obsp", par["obsp"]),
    ("varp", par["varp"]),
    ("uns", par["uns"]),
]

for slot_name, keys in _slots:
    if not keys:
        continue
    source_slot = getattr(source_mod, slot_name)
    target_slot = getattr(target_mod, slot_name)
    missing = [k for k in keys if k not in source_slot]
    if missing:
        raise ValueError(
            f"The following .{slot_name} keys were not found in source "
            f"modality '{par['source_modality']}': {missing}"
        )
    existing = [k for k in keys if k in target_slot]
    if existing and not par["allow_overwrite"]:
        raise ValueError(
            f"The following .{slot_name} keys already exist in the target "
            f"modality '{target_modality}': {existing}. "
            f"Use --allow_overwrite to overwrite them."
        )
    if existing:
        logger.warning("Overwriting existing .%s keys: %s", slot_name, existing)

    logger.info("Moving .%s keys: %s", slot_name, keys)
    for key in keys:
        target_slot[key] = source_slot[key]

logger.info(
    "Writing output to '%s' with compression '%s'",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    output_file=par["output"],
    h5mu=par["input_target"],
    modality_name=target_modality,
    modality_data=target_mod,
    output_compression=par["output_compression"],
)
