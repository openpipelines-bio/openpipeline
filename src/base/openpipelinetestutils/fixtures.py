from uuid import uuid4
import pytest

@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else f".{extension}"
        return tmp_path / f"{uuid4()}{extension}"
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