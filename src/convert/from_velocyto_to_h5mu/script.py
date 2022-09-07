from __future__ import annotations
from scanpy import read
from mudata import MuData
from pathlib import Path


## VIASH START
## VIASH END

def main():
    input_loom = par['input']
    if not Path(input_loom).is_file():
        raise FileNotFoundError(f"{input_loom} does not exist, is not a file or is not accessible.")
    adata = read(par['input'])
    result = MuData({'rna_velocity': adata})
    result.write(par['output'])
    

if __name__ == "__main__":
    main()