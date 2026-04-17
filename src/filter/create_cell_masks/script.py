import sys
from operator import eq, ge, gt, le, lt, ne

import mudata as md
import pandas as pd

################################################################################
# VIASH
################################################################################

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
}
## VIASH END

################################################################################
# FUNCTIONS
################################################################################


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


def parse_operator(operator_string):
    operators = {
        "lt": {"function": lt, "string": "lt", "symbol": "<"},
        "gt": {"function": gt, "string": "gt", "symbol": ">"},
        "le": {"function": le, "string": "le", "symbol": "<="},
        "ge": {"function": ge, "string": "ge", "symbol": ">="},
        "eq": {"function": eq, "string": "eq", "symbol": "=="},
        "ne": {"function": ne, "string": "ne", "symbol": "!="},
    }
    if operator_string not in operators:
        raise ValueError(
            "Operator must be one of 'lt', 'gt', 'le', 'ge', 'eq', or 'ne'. "
            f"Got: {operator_string}."
        )
    return operators[operator_string]


def parse_filters(raw_filters):
    if isinstance(raw_filters, str):
        raw_filters = [f for f in raw_filters.split(",") if f]

    filters = []
    for filter_string in raw_filters:
        parts = filter_string.split(":")
        if len(parts) not in {3, 4}:
            raise ValueError(
                "Each filter must be formatted as"
                "'<column>:<operator>:<value>:<group>' (<group> is optional)."
                f"Got: '{filter_string}'."
            )

        column, operator, value = parts[0], parts[1], parts[2]
        group = parts[3] if len(parts) == 4 else None

        operator = parse_operator(operator)

        filters.append(
            {
                "name": f"{group}_{column}_{operator['string']}_{value}"
                if group
                else f"{column}_{operator['string']}_{value}",
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
    masks = {}
    group_masks = {}
    overall_mask = pd.Series(True, index=adata.obs.index)

    for filter in filters:
        column = filter["column"]

        if column not in adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in adata.obs.")

        name = filter["name"]
        operator = filter["operator"]
        value = filter["value"]
        group = filter["group"]

        mask = operator(adata.obs[column], value)
        masks[name] = mask
        overall_mask &= mask

        if group:
            if group not in group_masks:
                group_masks[group] = pd.Series(True, index=adata.obs.index)
            group_masks[group] &= mask

    masks = pd.DataFrame(masks, index=adata.obs.index)
    group_masks = pd.DataFrame(group_masks, index=adata.obs.index)
    group_masks["overall"] = overall_mask

    return (masks, group_masks)


################################################################################
# MAIN
################################################################################


def main(par):
    print(f"====== Create cell masks (mudata v{md.__version__}) ======", flush=True)

    print(f"\n>>> Reading MuData from '{par['input']}'...", flush=True)
    mdata = md.read_h5mu(par["input"])
    print(mdata, flush=True)

    print(f"\n>>> Extracting modality '{par['modality']}'...", flush=True)
    if par["modality"] not in mdata.mod:
        raise KeyError(
            f"Modality '{par['modality']}' not found in MuData. "
            f"Available modalities: {list(mdata.mod.keys())}"
        )
    adata = mdata[par["modality"]]
    print(adata, flush=True)

    print("\n>>> Parsing filters...", flush=True)
    filters = parse_filters(par["filters"])
    print(f"Parsed {len(filters)} filters:", flush=True)
    for filter in filters:
        print(f"  - {filter['name']}: {filter['description']}", flush=True)

    print("\n>>> Creating masks...", flush=True)
    masks, group_masks = create_masks(adata, filters)
    print(f"Created {len(masks.columns)} individual masks", flush=True)
    print(masks, flush=True)
    print(f"\nCreated {len(group_masks.columns)} group masks", flush=True)
    print(group_masks, flush=True)

    print("\n>>> Adding masks to AnnData...", flush=True)
    obsm_name = f"{par['prefix']}_masks"
    adata.obsm[obsm_name] = masks
    print(f"Individual masks stored in obsm['{obsm_name}']", flush=True)
    print(adata.obsm[obsm_name], flush=True)

    group_mask_names = []
    for group in group_masks.columns:
        if group == "overall":
            mask_name = f"{par['prefix']}_mask"
        else:
            mask_name = f"{par['prefix']}_mask_{group}"

        adata.obs[mask_name] = group_masks[group]
        adata.obsm[obsm_name][group] = group_masks[group]
        group_mask_names.append(mask_name)

    print(f"\nGroup masks stored in obs with prefix '{par['prefix']}_mask'", flush=True)
    print(adata.obs[group_mask_names], flush=True)

    print("\n>>> Adding filters to AnnData...", flush=True)
    filters_name = f"{par['prefix']}_filters"
    filters_records = [
        {
            "name": filter["name"],
            "description": filter["description"],
            "column": filter["column"],
            "operator": filter["operator"].__name__,
            "value": filter["value"],
            "group": filter["group"],
        }
        for filter in filters
    ]
    adata.uns[filters_name] = pd.DataFrame(filters_records)
    print(f"Filters stored in uns['{filters_name}']", flush=True)
    print(adata.uns[filters_name], flush=True)

    print(f"\n>>> Writing output to '{par['output']}'...", flush=True)
    print(mdata, flush=True)
    mdata.write_h5mu(par["output"])

    print("\n>>> Done!\n")


if __name__ == "__main__":
    sys.exit(main(par))
