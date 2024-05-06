from mudata import read_h5mu
# import pytest
import sys
import pathlib

##VIASH START
par = {
    "input": "input.h5mu"
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

input_mudata = read_h5mu(par["input"])

def test_run():
    assert "rna" in input_mudata.mod.keys()
    assert input_mudata.n_obs == 713
    assert input_mudata.mod["rna"].var["feature_types"].unique() == [
        "Gene Expression"
    ], "Output X should only contain Gene Expression vars."
    
if __name__ == "__main__":
    with open("pytest.ini", "w") as ini_file:
        ini_file.write("[pytest]\npython_files=*.viash_script.py\nnorecursedirs=")
    
    import pytest
    sys.exit(pytest.main([str(pathlib.Path(__file__).resolve().parent),
                        #   "-o", "python_files=*.viash_script.sh",
                          "--verbose",
                        #   "--rootdir", str(pathlib.Path(__file__).resolve().parent)
                          ]))
    # sys.exit(pytest.main(["script.py"]))