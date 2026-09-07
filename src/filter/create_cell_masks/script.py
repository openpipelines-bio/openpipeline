import sys
from operator import eq, ge, gt, le, lt, ne

import pandas as pd
from mudata import read_h5ad

## VIASH START
par = {
    "input": "input.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "filters": [
        "total_counts:gt:1000:rna",
        "n_genes:gt:500:rna",
    ],
    "prefix": "cell",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402
from compress_h5mu import write_h5ad_to_h5mu_with_compression  # noqa: E402

logger = setup_logger()

OPERATORS = {
    "lt": {"function": lt, "string": "lt", "symbol": "<"},
    "gt": {"function": gt, "string": "gt", "symbol": ">"},
    "le": {"function": le, "string": "le", "symbol": "<="},
    "ge": {"function": ge, "string": "ge", "symbol": ">="},
    "eq": {"function": eq, "string": "eq", "symbol": "=="},
    "ne": {"function": ne, "string": "ne", "symbol": "!="},
}


def parse_value(raw_value):
    if raw_value.lower() in {"true", "false"}:
        return raw_value.lower() == "true"
    try:
        return int(raw_value)
    except ValueError:
        pass
    try:
        return float(raw_value)
    except ValueError:
        pass
    return raw_value


def parse_operator(operator_string, filter_string):
    if operator_string not in OPERATORS:
        raise ValueError(
            f"Unknown operator '{operator_string}' in filter '{filter_string}'. "
            f"Must be one of: {', '.join(OPERATORS)}."
        )
    return OPERATORS[operator_string]


def parse_filters(raw_filters):
    filters = []
    for filter_string in raw_filters:
        parts = filter_string.split(":")
        if len(parts) not in {3, 4}:
            raise ValueError(
                f"Each filter must be formatted as "
                f"'<column>:<operator>:<value>:<group>' (<group> is optional). "
                f"Got: '{filter_string}'."
            )

        column, operator_str, value = parts[0], parts[1], parts[2]
        group = parts[3] if len(parts) == 4 else None
        operator = parse_operator(operator_str, filter_string)

        name_parts = [p for p in (group, column, operator["string"], value) if p]
        filters.append(
            {
                "name": "_".join(name_parts),
                "description": f"{column} {operator['symbol']} {value}"
                + (f" ({group})" if group else ""),
                "column": column,
                "operator": operator["function"],
                "value": parse_value(value),
                "group": group,
            }
        )

    return filters


def create_masks(adata, filters):
    missing = sorted({f["column"] for f in filters} - set(adata.obs.columns))
    if missing:
        raise KeyError(
            f"The following columns referenced by filters are not in .obs: {missing}"
        )

    masks = {}
    group_masks = {}
    overall_mask = pd.Series(True, index=adata.obs.index)

    for filt in filters:
        mask = filt["operator"](adata.obs[filt["column"]], filt["value"])
        masks[filt["name"]] = mask
        overall_mask &= mask

        group = filt["group"]
        if group:
            if group not in group_masks:
                group_masks[group] = pd.Series(True, index=adata.obs.index)
            group_masks[group] &= mask

    masks = pd.DataFrame(masks, index=adata.obs.index)
    group_masks = pd.DataFrame(group_masks, index=adata.obs.index)
    group_masks["overall"] = overall_mask

    return masks, group_masks


def main(par):
    prefix = par["prefix"] or ""
    prefix_part = f"{prefix}_" if prefix else ""

    logger.info("Reading modality '%s' from '%s'", par["modality"], par["input"])
    try:
        adata = read_h5ad(par["input"], mod=par["modality"])
    except KeyError:
        raise ValueError(f"Modality '{par['modality']}' not found in '{par['input']}'.")

    logger.info("Parsing %d filter(s)", len(par["filters"]))
    filters = parse_filters(par["filters"])
    for filt in filters:
        logger.info("  - %s: %s", filt["name"], filt["description"])

    logger.info("Creating masks")
    masks, group_masks = create_masks(adata, filters)
    logger.info(
        "Created %d individual mask(s) and %d group mask(s)",
        len(masks.columns),
        len(group_masks.columns),
    )

    obsm_name = f"{prefix_part}masks"
    adata.obsm[obsm_name] = masks
    logger.info("Stored individual masks in .obsm['%s']", obsm_name)

    group_mask_names = []
    for group in group_masks.columns:
        mask_suffix = "" if group == "overall" else f"_{group}"
        mask_name = f"{prefix_part}mask{mask_suffix}"
        adata.obs[mask_name] = group_masks[group]
        adata.obsm[obsm_name][group] = group_masks[group]
        group_mask_names.append(mask_name)
    logger.info("Stored group masks in .obs: %s", group_mask_names)

    filters_name = f"{prefix_part}filters"
    filters_records = [
        {
            "name": filt["name"],
            "description": filt["description"],
            "column": filt["column"],
            "operator": filt["operator"].__name__,
            "value": filt["value"],
            "group": filt["group"],
        }
        for filt in filters
    ]
    filters_df = pd.DataFrame(filters_records)
    # Empty string for ungrouped filters: h5 cannot write Python None as a
    # string and anndata does not opt into nullable string writing by default.
    filters_df["group"] = filters_df["group"].fillna("").astype(str)
    adata.uns[filters_name] = filters_df
    logger.info("Stored filter definitions in .uns['%s']", filters_name)

    logger.info(
        "Writing output to '%s' with compression '%s'",
        par["output"],
        par["output_compression"],
    )
    write_h5ad_to_h5mu_with_compression(
        output_file=par["output"],
        h5mu=par["input"],
        modality_name=par["modality"],
        modality_data=adata,
        output_compression=par["output_compression"],
    )


if __name__ == "__main__":
    sys.exit(main(par))
