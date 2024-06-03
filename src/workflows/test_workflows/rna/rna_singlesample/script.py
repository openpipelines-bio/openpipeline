from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

##VIASH START
par = {
    "input": "input.h5mu",
    "input_og": "input_og.h5mu",
    "input_id": "mitochondrial_test"
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

input_mudata = read_h5mu(par["input_og"])
output_mudata = read_h5mu(par["input"])

# General assertions for all cases
assert input_mudata.mod.keys() == output_mudata.mod.keys(), "Input and output should have the same modalities."
assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])
assert set(['filter_with_counts',
            'scrublet_doublet_score',
            'filter_with_scrublet']).issubset(output_mudata.mod["rna"].obs)
assert 'filter_with_counts' in output_mudata.mod["rna"].var

# Specific assertions for each test case
match par["input_id"]:
    case "mitochondrial_test":
        assert output_mudata.mod["rna"].n_obs < input_mudata.mod["rna"].n_obs, "Number of observations should be less than the input."
        assert set(['fraction_mitochondrial',
                    'total_counts_mitochondrial',
                    'pct_mitochondrial',
                    'filter_mitochondrial']).issubset(output_mudata.mod["rna"].obs.keys())
        assert 'mitochondrial' in output_mudata.mod["rna"].var.keys()
        pass
    case "simple_execution_test":
        pass
    case "test_different_fraction_column":
        assert output_mudata.mod["rna"].n_obs < input_mudata.mod["rna"].n_obs, "Number of observations should be less than the input."
        assert set(['foobar',
                    'total_counts_mitochondrial',
                    'pct_mitochondrial',
                    'filter_mitochondrial']).issubset(output_mudata.mod["rna"].obs.keys()) and\
                        'fraction_mitochondrial' not in output_mudata.mod["rna"].obs.keys()
        assert 'mitochondrial' in output_mudata.mod["rna"].var.keys()
        pass