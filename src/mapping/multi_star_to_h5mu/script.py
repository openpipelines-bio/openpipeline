from pathlib import Path
import pandas as pd
import mudata as md
import anndata as ad
import numpy as np
import json

## VIASH START
par = {"input": "output/A2_raw", "output": "test_output.h5mu"}
meta = {"temp_dir": "/tmp"}
## VIASH END

# convert to path
input_dir = Path(par["input"])

# read counts information
print("> Read counts data", flush=True)
per_obs_data = []

for input_counts in (input_dir / "per").glob("**/htseq-count.txt"):
    per_obs_dir = input_counts.parent
    input_id = per_obs_dir.name
    input_multiqc = per_obs_dir / "multiqc_data" / "multiqc_data.json"

    data = pd.read_table(
        input_counts,
        index_col=0,
        names=["cell_id", input_id],
        dtype={"cell_id": "U", input_id: "i"},
    )
    data2 = data[~data.index.str.startswith("__")]

    with open(input_multiqc, "r") as file:
        qc = json.load(file)

    qc_star = qc.get("report_saved_raw_data", {}).get("multiqc_star", {}).get(input_id)
    qc_htseq = (
        qc.get("report_saved_raw_data", {}).get("multiqc_htseq", {}).get("htseq-count")
    )

    per_obs_data.append(
        {
            "counts": data2.transpose(),
            "qc_star": pd.DataFrame(qc_star, index=[input_id]),
            "qc_htseq": pd.DataFrame(qc_htseq, index=[input_id]),
        }
    )


# combine all counts
counts = pd.concat([x["counts"] for x in per_obs_data], axis=0)
qc_star = pd.concat([x["qc_star"] for x in per_obs_data], axis=0)
qc_htseq = pd.concat([x["qc_htseq"] for x in per_obs_data], axis=0)

# read feature info
feature_info = pd.read_csv(input_dir / "feature_info.tsv", sep="\t", index_col=0)
feature_info_ord = feature_info.loc[counts.columns]

var = pd.DataFrame(
    data={
        "gene_ids": feature_info_ord.index,
        "feature_types": "Gene Expression",
        "gene_name": feature_info_ord["feature_name"],
    }
).set_index("gene_ids")

print("> construct anndata", flush=True)
adata = ad.AnnData(
    X=counts, obsm={"qc_star": qc_star, "qc_htseq": qc_htseq}, var=var, dtype=np.int32
)

print("> convert to mudata", flush=True)
mdata = md.MuData(adata)

print("> write to file", flush=True)
mdata.write_h5mu(par["output"], compression=par["output_compression"])
