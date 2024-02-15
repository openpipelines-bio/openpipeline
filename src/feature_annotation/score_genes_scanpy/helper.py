from typing import List, Dict, Any, Optional

def read_gene_list(
        par: Dict[str, Any],
        gene_names: List[str],
        list_key: str,
        file_key: str,
        required: bool = True) -> Optional[List[str]]:
    """
    Reads a gene list from the parameters and returns it as a list of strings.
    """

    # check whether one or the other was provided, if required
    if required and not par[list_key] and not par[file_key]:
        raise ValueError(f"Either --{list_key} or --{file_key} must be set")

    # read gene list from parameters
    list = par[list_key] if par[list_key] else []

    # read gene list from file
    if par[file_key]:
        with open(par[file_key]) as file:
            file_genes = [x.strip() for x in file]
        list.extend(file_genes)

    # check for missing genes
    if not par["allow_missing_genes"] and list:
        missing = set(list).difference(gene_names)
        if missing:
            raise ValueError(f"The follow genes are missing from the input dataset: {missing}")

    # return gene list
    if required:
        if not list:
            raise ValueError(f"No genes detected in --{list_key} or --{file_key}")
        return list
    else:
        # turn empty list into None
        if not list:
            return None
        return list
