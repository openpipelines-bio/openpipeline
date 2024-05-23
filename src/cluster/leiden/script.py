import sys
import signal
import os
import threading
import time
import mudata as mu
import pandas as pd
import scanpy as sc
import numpy as np
import numpy.typing as npt
import anndata as ad
from multiprocessing import managers, shared_memory, get_context
from concurrent.futures import ProcessPoolExecutor, process
from scipy.sparse import csr_matrix
from pathlib import Path
from itertools import repeat
import shutil

## VIASH START
par = {
    "input": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "output_format": "h5mu",
    "obsm_name": "leiden",
    "resolution": [1, 0.25, 0.10, 0.05, 0.01, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.7, 0.3, 0.35, 0.95],
    "obsp_connectivities": "connectivities",
    "uns_name": "leiden",
    "output_compression": "gzip"
}
meta = {
    "cpus": 10,
    "resources_dir": '.'
}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger

# START TEMPORARY WORKAROUND compress_h5mu
# reason: resources aren't available when using Nextflow fusion
# from compress_h5mu import compress_h5mu
from h5py import File as H5File
from h5py import Group, Dataset
from typing import Union, Literal
from functools import partial

def compress_h5mu(input_path: Union[str, Path], 
                output_path: Union[str, Path], 
                compression: Union[Literal['gzip'], Literal['lzf']]):
    input_path, output_path = str(input_path), str(output_path)

    def copy_attributes(in_object, out_object):
        for key, value in in_object.attrs.items():
            out_object.attrs[key] = value

    def visit_path(output_h5: H5File,
                   compression: Union[Literal['gzip'], Literal['lzf']], 
                   name: str, object: Union[Group, Dataset]):
            if isinstance(object, Group):
                new_group = output_h5.create_group(name)
                copy_attributes(object, new_group)
            elif isinstance(object, Dataset):
                # Compression only works for non-scalar Dataset objects
                # Scalar objects dont have a shape defined
                if not object.compression and object.shape not in [None, ()]: 
                    new_dataset = output_h5.create_dataset(name, data=object, compression=compression)
                    copy_attributes(object, new_dataset)
                else:
                    output_h5.copy(object, name)
            else:
                raise NotImplementedError(f"Could not copy element {name}, "
                                          f"type has not been implemented yet: {type(object)}")

    with H5File(input_path, 'r') as input_h5, H5File(output_path, 'w', userblock_size=512) as output_h5:
        copy_attributes(input_h5, output_h5)
        input_h5.visititems(partial(visit_path, output_h5, compression))

    with open(input_path, "rb") as input_bytes:
        # Mudata puts metadata like this in the first 512 bytes:
        # MuData (format-version=0.1.0;creator=muon;creator-version=0.2.0)
        # See mudata/_core/io.py, read_h5mu() function
        starting_metadata = input_bytes.read(100)
        # The metadata is padded with extra null bytes up until 512 bytes
        truncate_location = starting_metadata.find(b"\x00")
        starting_metadata = starting_metadata[:truncate_location]
    with open(output_path, "br+") as f:
        nbytes = f.write(starting_metadata)
        f.write(b"\0" * (512 - nbytes)) 
# END TEMPORARY WORKAROUND compress_h5mu

logger = setup_logger()

class SharedNumpyMatrix():
    def __init__(self, shared_memory: shared_memory.SharedMemory, dtype: npt.DTypeLike, shape: tuple[int, int]) -> None:
        self._memory = shared_memory
        self._dtype = dtype
        self._shape = shape
    
    @classmethod
    def from_numpy(cls, memory_manager: managers.SharedMemoryManager, array: npt.ArrayLike):
        shm = memory_manager.SharedMemory(size=array.nbytes)
        array_in_shared_memory = np.ndarray(array.shape, dtype=array.dtype, buffer=shm.buf)
        # Copy the data into shared memory
        array_in_shared_memory[:] = array[:]
        return cls(shm, array.dtype, array.shape)

    def to_numpy(self):
       return np.ndarray(self._shape, dtype=self._dtype, buffer=self._memory.buf)
    
    def close(self):
        self._memory.close()

class SharedCsrMatrix():
    def __init__(self,
                 data: SharedNumpyMatrix, 
                 indices: SharedNumpyMatrix, 
                 indptr: SharedNumpyMatrix, 
                 shape: npt.DTypeLike):
        self._data = data
        self._indices = indices
        self._indptr = indptr
        self._shape = shape

    @classmethod
    def from_csr_matrix(cls, memory_manager: managers.SharedMemoryManager, csr_matrix_obj: csr_matrix):
        return cls(
            SharedNumpyMatrix.from_numpy(memory_manager, csr_matrix_obj.data),
            SharedNumpyMatrix.from_numpy(memory_manager, csr_matrix_obj.indices),
            SharedNumpyMatrix.from_numpy(memory_manager, csr_matrix_obj.indptr),
            csr_matrix_obj.shape,
        )
    
    def to_csr_matrix(self):
        return csr_matrix(
            (self._data.to_numpy(), self._indices.to_numpy(), self._indptr.to_numpy()),
            shape=self._shape, 
            copy=False
        )

    def close(self):
        self._data.close()
        self._indices.close()
        self._indptr.close()

def create_empty_anndata_with_connectivities(connectivities, obs_names):
    empty_anndata = ad.AnnData(np.zeros((connectivities.shape[0], 1)),
                               obs=pd.DataFrame(index=list(obs_names)))
    empty_anndata.obsp['connectivities'] = connectivities
    return empty_anndata

def run_single_resolution(shared_csr_matrix, obs_names, resolution):
    try:
        connectivities = shared_csr_matrix.to_csr_matrix()
        adata = create_empty_anndata_with_connectivities(connectivities, obs_names)
        adata_out = sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=str(resolution),
            obsp="connectivities",
            copy=True
            )
        return adata_out.obs[str(resolution)]
    finally:
        obs_names.shm.close()
        shared_csr_matrix.close() 

def start_child_process_terminator(ppid):
    pid = os.getpid()

    def terminator():
        while True:
            try:
                # If sig is 0, then no signal is sent, but error checking is still performed; 
                # this can be used to check for the existence of a process ID
                os.kill(ppid, 0)
            except OSError:
                os.kill(pid, signal.SIGTERM)
            time.sleep(1)

    thread = threading.Thread(target=terminator, daemon=True)
    thread.start()

def main():
    logger.info("Reading %s.", par["input"])
    adata = mu.read_h5ad(par["input"], mod=par['modality'], backed='r')
    logger.info("Processing modality '%s'.", par['modality'])
    try:
        connectivities = adata.obsp[par['obsp_connectivities']]
    except KeyError:
        raise ValueError(f"Could not find .obsp key \"{par['obsp_connectivities']}\" "
                        "in modality {par['modality']}")

    with managers.SharedMemoryManager() as smm:
        # anndata converts the index to strings, so no worries that it cannot be stored in ShareableList
        # because it has an unsupported dtype. It should always be string...
        index_contents = adata.obs.index.to_list()
        assert all([isinstance(item, str) for item in index_contents])
        obs_names = smm.ShareableList(index_contents)

        shared_csr_matrix = SharedCsrMatrix.from_csr_matrix(smm, connectivities)
        exit_with_other_code = None
        with ProcessPoolExecutor(max_workers=meta['cpus'], max_tasks_per_child=1, 
                                 mp_context=get_context('spawn'),
                                 initializer=start_child_process_terminator,
                                 initargs=(os.getpid(),)) as executor:
            results = executor.map(run_single_resolution, 
                                    repeat(shared_csr_matrix), 
                                    repeat(obs_names), 
                                    par["resolution"],
                                    chunksize=1)
            try:
                results = {str(resolution): result for resolution, result 
                        in zip(par["resolution"], results)}
            except process.BrokenProcessPool as e:
                # This assumes that one of the child processses was killed by the kernel
                # because the oom killer was activated. This the is the most likely scenario,
                # other causes could be:
                # * Subprocess terminates without raising a proper exception.
                # * The code of the process handling the communication is broke (i.e. a python bug)
                # * The return data could not be pickled.
                print(e, file=sys.stderr, flush=True)
                exit_with_other_code = 137
                raise e
        shared_csr_matrix.close()
        obs_names.close() 
    if exit_with_other_code:
        from multiprocessing import resource_tracker
        resource_tracker._resource_tracker._stop()
        sys.exit(exit_with_other_code)

    adata.obsm[par["obsm_name"]] = pd.DataFrame(results)
    logger.info("Writing to %s.", par["output"])

    output_file = Path(par["output"])
    logger.info('Writing output to %s.', par['output'])
    output_file_uncompressed = output_file.with_name(output_file.stem + "_uncompressed.h5mu") \
        if par["output_compression"] else output_file
    shutil.copyfile(par['input'], output_file_uncompressed)
    mu.write_h5ad(filename=output_file_uncompressed, mod=par['modality'], data=adata)
    if par["output_compression"]:
        compress_h5mu(output_file_uncompressed, output_file, compression=par["output_compression"])
        output_file_uncompressed.unlink()
    logger.info("Finished.")

if __name__ == "__main__":
    main()