from muon import prot as pt
from mudata import read_h5mu

## VIASH START
par = {
    'input': '/home/di/code/openpipeline/work/f1/3a5269aaf10aaf09fee85b0d4bfcab/pbmc.do_filter.output.h5mu',
    'modality': 'prot',
    'output': "foo.h5mu"
}
## VIASH END


def main():
    input_h5mu = read_h5mu(par['input'])
    modality = input_h5mu[par['modality']]
    normalized_counts = pt.pp.clr(modality, inplace=False)
    input_h5mu[par["modality"]].layers['clr'] = normalized_counts.X
    input_h5mu.write_h5mu(par['output'], compression=par["output_compression"])

if __name__ == "__main__":
    main()