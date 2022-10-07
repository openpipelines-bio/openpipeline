import mudata
import popv
from torch.cuda import is_available as cuda_is_available
try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
    }
### VIASH END

