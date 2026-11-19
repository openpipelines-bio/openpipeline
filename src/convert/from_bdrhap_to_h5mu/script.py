import mudata as mu
import pandas as pd
import scirpy as ir

## VIASH START
par = {
    "id": "sample",
    "input": "resources_test/bdrhap_5kjrt/processed/output_raw/sample.cellismo",
    "input_vdj_contigs": None,
    "output": "bd_rhap_to_h5mu_test.h5mu",
    "output_compression": None,
}
## VIASH END

# Modality names for the TCR/BCR loci, consistent with the VDJ modalities
# created by convert/from_cellranger_multi_to_h5mu
LOCUS_TO_MODALITY = {
    "TRA": "vdj_t",
    "TRB": "vdj_t",
    "TRG": "vdj_t_gd",
    "TRD": "vdj_t_gd",
    "IGH": "vdj_b",
    "IGK": "vdj_b",
    "IGL": "vdj_b",
}

print(">> Reading input file", flush=True)
mdata = mu.read_h5mu(par["input"])

# Check if modalities are present
modalities = list(mdata.mod.keys())
assert len(modalities) > 0, "No modalities found in input data"

# BD Rhapsody >= v4.0 stores the dominant TCR/BCR contigs of the putative cells in
# the "airr" modality. The VDJ modalities are created from the unfiltered contigs
# instead, which contain all contigs for all cell labels.
if "airr" in mdata.mod and not par["input_vdj_contigs"]:
    raise ValueError(
        "The input contains VDJ data (modality 'airr'). Please also provide the "
        "unfiltered VDJ contigs ('*_VDJ_Unfiltered_Contigs_AIRR.tsv.gz') using "
        "--input_vdj_contigs."
    )


def process_modality_inline(adata, modality, pipeline_inputs):
    adata.obs["library_id"] = " & ".join(pipeline_inputs["Libraries"])
    adata.obs["cell_id"] = adata.obs.index
    adata.obs["run_id"] = par["id"]

    # 'sample_id' is not used for the sample name, because it is overwritten
    # with the sample id of the input when processing multiple samples.
    adata.obs.rename(
        columns={"Sample_Tag": "sample_tag", "Sample_Name": "sample_name"},
        inplace=True,
    )

    # The dimensionality reductions calculated by the BD Rhapsody pipeline are
    # renamed to avoid them being confused with (or overwritten by) dimensionality
    # reductions calculated downstream.
    for key in ["X_tsne", "X_umap"]:
        if key in adata.obsm:
            adata.obsm[f"{key}_bd_rhapsody"] = adata.obsm.pop(key)

    # The VDJ modalities do not contain features
    if modality in LOCUS_TO_MODALITY.values():
        return

    # Use the gene ids as index and store the gene symbols in .var['gene_symbol'].
    # BD Rhapsody >= v3.0 stores the gene ids (e.g. Ensembl ids) in .var["gene_id"].
    # For targeted assays and AbSeq, the gene id is equal to the feature name.
    gene_symbols = (
        adata.var.pop("gene_name") if "gene_name" in adata.var else adata.var.index
    )
    adata.var["gene_symbol"] = gene_symbols.to_numpy()
    if "gene_id" in adata.var:
        gene_ids = adata.var.pop("gene_id")
        if not gene_ids.is_unique:
            raise ValueError(
                f"The gene ids of modality '{modality}' are not unique: "
                f"{gene_ids[gene_ids.duplicated()].unique().tolist()[:10]}"
            )
        adata.var.index = gene_ids.to_numpy()
    adata.var.index.name = None

    if modality == "rna":
        adata.var["feature_types"] = "Gene Expression"
        # WTA runs use a reference archive, targeted runs use targeted reference fasta files
        if "Reference_Archive" in pipeline_inputs:
            adata.var["reference_file"] = pipeline_inputs["Reference_Archive"]
        else:
            adata.var["reference_file"] = " & ".join(
                pipeline_inputs["Targeted_Reference"]
            )

    elif modality == "prot":
        adata.var["feature_types"] = "Antibody Capture"
        adata.var["reference_file"] = " & ".join(pipeline_inputs["AbSeq_Reference"])

    # TODO: add other modalities


def process_vdj_contigs(contigs_file):
    # Create one modality per receptor type, containing all contigs of all cell labels.
    # The 'putative_cell' column indicates whether a cell label is a putative cell.
    contigs = pd.read_csv(
        contigs_file, sep="\t", dtype={"cell_id": str}, low_memory=False
    )
    unknown_loci = set(contigs["locus"]) - set(LOCUS_TO_MODALITY)
    assert not unknown_loci, f"Unexpected loci in VDJ contigs: {unknown_loci}"

    vdj_modalities = {}
    for modality, contigs_modality in contigs.groupby(
        contigs["locus"].map(LOCUS_TO_MODALITY)
    ):
        print(f">> Creating modality '{modality}' from VDJ contigs", flush=True)
        adata = ir.io.read_airr(contigs_modality, cell_attributes=["putative_cell"])
        # scirpy names the index 'cell_id', which clashes with the 'cell_id' column
        # once the index is altered (e.g. prefixed with the sample id by add_id).
        adata.obs.index.name = None
        vdj_modalities[modality] = adata
    return vdj_modalities


# BD Rhapsody >= v2.3 stores the pipeline inputs in the global .uns,
# older versions store them in the .uns of each modality
pipeline_inputs = mdata.uns.get("Pipeline_Inputs") or next(
    (
        modality.uns["Pipeline_Inputs"]
        for modality in mdata.mod.values()
        if "Pipeline_Inputs" in modality.uns
    ),
    None,
)
assert pipeline_inputs is not None, "Pipeline_Inputs not found in input data"

result_modalities = {key: value for key, value in mdata.mod.items() if key != "airr"}
if par["input_vdj_contigs"]:
    result_modalities |= process_vdj_contigs(par["input_vdj_contigs"])

for key, value in result_modalities.items():
    print(">> Processing modality:", key, flush=True)
    process_modality_inline(value, key, pipeline_inputs)

output = mu.MuData(result_modalities, uns=mdata.uns)

print(">> Writing output file", flush=True)
output.write_h5mu(par["output"], compression=par["output_compression"])
