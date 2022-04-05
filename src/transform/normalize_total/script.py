import scanpy as sc
import muon as mu

## VIASH START
par = {
    "input": "work/d9/3adbd080e0de618d44b59b1ec81685/run.output.h5mu",
    "output": "output.h5mu",
    "target_sum": 10000,
    "modality": ["rna"],
    "exclude_highly_expressed": False
}
meta = {"functionality_name": "lognorm"}
## VIASH END

print("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

print(par)

for mod in par["modality"]:
    print(f"Performing total normalization on modality {mod}")
    dat = mdata.mod[mod]
    print(dat)
    sc.pp.normalize_total(
        dat,
        # target_sum=par["target_sum"],
        # exclude_highly_expressed=par["exclude_highly_expressed"]
    )

print("Writing to file")
mdata.write(filename=par["output"])
