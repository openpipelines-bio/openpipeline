import zarr
import h5py
from contextlib import contextmanager


@contextmanager
def mudata_opener(file_loc, mode=None):
    open_mudata = None
    input_is_zarr = False
    try:
        open_mudata = zarr.open(file_loc, zarr_format=3, mode=mode)
        input_is_zarr = True
        yield open_mudata, input_is_zarr
    except (zarr.errors.GroupNotFoundError, NotADirectoryError):
        try:
            open_mudata = h5py.File(file_loc, mode=mode)
            yield open_mudata, input_is_zarr
        except (FileNotFoundError, IsADirectoryError, KeyError) as e:
            e.add_note(f"Could not open file {file_loc}.")
            raise e
        finally:
            try:
                if open_mudata:
                    open_mudata.close()
                    del open_mudata
            except (AttributeError, UnboundLocalError):
                pass
