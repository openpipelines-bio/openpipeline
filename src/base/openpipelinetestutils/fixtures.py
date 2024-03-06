from uuid import uuid4
import pytest
import pandas as pd
import anndata as ad
import mudata as md

@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else extension
        return tmp_path / f"{uuid4()}.{extension}"
    return wrapper 

@pytest.fixture
def random_h5mu_path(random_path):
    def wrapper():
        return random_path(extension="h5mu")
    return wrapper

@pytest.fixture
def write_mudata_to_file(random_h5mu_path):
    def wrapper(mudata_obj):
        output_path = random_h5mu_path()
        mudata_obj.write(output_path)
        return output_path
    return wrapper

@pytest.fixture
def small_anndata_1():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]],
                        index=df.columns, columns=["Feat"])
    ad1 = ad.AnnData(df, obs=obs, var=var)
    return ad1

@pytest.fixture
def small_anndata_2():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = ad.AnnData(df, obs=obs2, var=var2)
    return ad2

@pytest.fixture
def small_mudata(small_anndata_1, small_anndata_2):
    return md.MuData({'mod1': small_anndata_1, 'mod2': small_anndata_2})

@pytest.fixture
def small_mudata_path(small_mudata, write_mudata_to_file):
    return write_mudata_to_file(small_mudata)

@pytest.fixture
def split_small_mudata_path(small_mudata_mod1_path, small_mudata_mod2_path):
    return small_mudata_mod1_path, small_mudata_mod2_path

@pytest.fixture
def small_mudata_mod1_path(small_mudata, write_mudata_to_file):
    return write_mudata_to_file(md.MuData({'mod1': small_mudata.mod['mod1']}))

@pytest.fixture
def small_mudata_mod2_path(small_mudata, write_mudata_to_file):
    return write_mudata_to_file(md.MuData({'mod2': small_mudata.mod['mod2']}))

