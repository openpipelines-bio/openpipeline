import sys

import mudata

from cellphonedb.utils.search_utils import search_analysis_results

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "output_compression": "gzip",
    "modality": "rna",
    "uns_significant_means": "cellphonedb_significant_means",
    "uns_deconvoluted": "cellphonedb_deconvoluted",
    "uns_interaction_scores": "cellphonedb_interaction_scores",
    "query_cell_types_1": None,
    "query_cell_types_2": None,
    "query_genes": None,
    "query_interactions": None,
    "query_classifications": None,
    "query_minimum_score": None,
    "separator": "|",
    "long_format": False,
    "output_uns_search_results": "cellphonedb_search_results",
}
meta = {"resources_dir": "."}
### VIASH END

sys.path.append(meta["resources_dir"])
from compress_h5mu import write_h5ad_to_h5mu_with_compression


def sanitize_dataframe_for_h5ad(df):
    # CellphoneDB leaves NaN in object-dtype columns such as gene_a/gene_b for
    # complex partners (which have no single gene name); h5ad can't write a
    # column that mixes strings and floats.
    df = df.copy()
    object_columns = df.select_dtypes(include="object").columns
    df[object_columns] = df[object_columns].fillna("").astype(str)
    return df


def main():
    mod = mudata.read_h5ad(par["input"], mod=par["modality"])

    if par["uns_significant_means"] not in mod.uns:
        raise ValueError(
            f"'{par['uns_significant_means']}' not found in .uns for modality "
            f"{par['modality']}. Run interpret/cellphonedb first."
        )
    if par["uns_deconvoluted"] not in mod.uns:
        raise ValueError(
            f"'{par['uns_deconvoluted']}' not found in .uns for modality "
            f"{par['modality']}. Run interpret/cellphonedb first."
        )

    interaction_scores = None
    if par["query_minimum_score"] is not None:
        if par["uns_interaction_scores"] not in mod.uns:
            raise ValueError(
                "--query_minimum_score requires interaction scores in .uns "
                f"(key '{par['uns_interaction_scores']}'); rerun interpret/cellphonedb "
                "with --score_interactions to compute them."
            )
        interaction_scores = mod.uns[par["uns_interaction_scores"]]

    result = search_analysis_results(
        query_cell_types_1=par["query_cell_types_1"],
        query_cell_types_2=par["query_cell_types_2"],
        query_genes=par["query_genes"],
        query_interactions=par["query_interactions"],
        query_classifications=par["query_classifications"],
        query_minimum_score=par["query_minimum_score"],
        significant_means=mod.uns[par["uns_significant_means"]],
        deconvoluted=mod.uns[par["uns_deconvoluted"]],
        interaction_scores=interaction_scores,
        separator=par["separator"],
        long_format=par["long_format"],
    )

    if result is None:
        raise ValueError(
            "CellphoneDB search_analysis_results returned no result; check that "
            f"'{par['uns_significant_means']}' and '{par['uns_deconvoluted']}' contain "
            "valid CellphoneDB output tables."
        )

    mod.uns[par["output_uns_search_results"]] = sanitize_dataframe_for_h5ad(result)

    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], mod, par["output_compression"]
    )


if __name__ == "__main__":
    main()
