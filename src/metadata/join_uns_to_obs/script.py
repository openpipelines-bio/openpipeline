from mudata import read_h5mu
import pandas as pd

import logging
from sys import stdout

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

### VIASH START
par = {
    "input": "work/f5/5f6365898ca5a42a360301a0c9e200/TSP15_Eye_ScleraEtc_10X_2_1.add_id.output.h5mu",
    "uns_key": "data",
    "output": "foo.h5mu",
    "modality": "rna"
}
### VIASH END


logger.info("Read mudata from file")
mdata = read_h5mu(par['input'])
mod_data = mdata.mod[par['modality']]

logger.info("Joining uns to obs")
mod_data.obs = mod_data.obs.join(mod_data.uns[par['uns_key']])

logger.info("Write output to mudata file")
mdata.write_h5mu(par['output'])

        

