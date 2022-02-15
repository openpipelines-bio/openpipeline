import scanpy as sc
import muon as mu
import multiprocessing

### VIASH START
par = {
  'input': 'resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'output': 'output.h5mu',
  'normalized_umi_count': 10000,
  'regress_out_variables': [ ],
  'modality': [ 'rna' ]
}
meta = { 'functionality_name': 'lognorm' }
### VIASH END

print("Reading input mudata")
mudata = mu.read_h5mu(par["input"])

for mod in par['modality']:
    print(f"Performing log normalization on modality {mod}")
    data = mudata.mod[mod]
    
    sc.pp.normalize_total(data, target_sum = par["normalized_umi_count"])
    sc.pp.log1p(data)

    if par["regress_out_variables"] is not None and len(par["regress_out_variables"]) > 0:
        print("Regress out variables on modality {mod}")        
        sc.pp.regress_out(data, par["regress_out_variables"], n_jobs=multiprocessing.cpu_count()-1)

# can we assume execution_log exists?
if mudata.uns is None or "execution_log" not in mudata.uns:
    mudata.uns["execution_log"] = []
# store new entry
new_entry = { "component": meta["functionality_name"], "params": par }
mudata.uns["execution_log"] = mudata.uns["execution_log"] + [ new_entry ]

print("Writing to file")
mudata.write(par["output"])