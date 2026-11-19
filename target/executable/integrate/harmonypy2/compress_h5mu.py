import shutil
from anndata import AnnData
from mudata import write_h5ad
from h5py import File as H5File
from h5py import Group, Dataset
from pathlib import Path
from typing import Union, Literal
from functools import partial


def compress_h5mu(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    compression: Union[Literal["gzip"], Literal["lzf"]],
):
    input_path, output_path = str(input_path), str(output_path)

    def copy_attributes(in_object, out_object):
        for key, value in in_object.attrs.items():
            out_object.attrs[key] = value

    def visit_path(
        output_h5: H5File,
        compression: Union[Literal["gzip"], Literal["lzf"]],
        name: str,
        object: Union[Group, Dataset],
    ):
        if isinstance(object, Group):
            new_group = output_h5.create_group(name)
            copy_attributes(object, new_group)
        elif isinstance(object, Dataset):
            # Compression only works for non-scalar Dataset objects
            # Scalar objects dont have a shape defined
            if not object.compression and object.shape not in [None, ()]:
                new_dataset = output_h5.create_dataset(
                    name, data=object, compression=compression
                )
                copy_attributes(object, new_dataset)
            else:
                output_h5.copy(object, name)
        else:
            raise NotImplementedError(
                f"Could not copy element {name}, "
                f"type has not been implemented yet: {type(object)}"
            )

    with (
        H5File(input_path, "r") as input_h5,
        H5File(output_path, "w", userblock_size=512) as output_h5,
    ):
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


def write_h5ad_to_h5mu_with_compression(
    output_file: Union[str, Path],
    h5mu: Union[str, Path],
    modality_name: str,
    modality_data: AnnData,
    output_compression=None,
):
    output_file = Path(output_file)
    h5mu = Path(h5mu)
    output_file_uncompressed = (
        output_file.with_name(output_file.stem + "_uncompressed.h5mu")
        if output_compression
        else output_file
    )
    shutil.copyfile(h5mu, output_file_uncompressed)
    write_h5ad(filename=output_file_uncompressed, mod=modality_name, data=modality_data)
    if output_compression:
        compress_h5mu(
            output_file_uncompressed, output_file, compression=output_compression
        )
        output_file_uncompressed.unlink()
