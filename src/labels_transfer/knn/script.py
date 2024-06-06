import mudata as mu
import pandas as pd
import re
from pathlib import Path

### VIASH START
par = {
  'input': 'reference_download/reference.h5mu',
  'modality': 'rna',
  'obs_feature': 'donor_assay',
  'output': 'reference_download/sample_split',
  'output_compression': None,
  'output_files': 'reference_download/sample_split/sample_files.csv'
}

### VIASH END

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
logger = setup_logger()


def main():
    logger.info(f"Reading {par['input']}")
    input_file = Path(par["input"].strip())

    mdata = mu.read_h5mu(input_file)
    adata = mdata.mod[par["modality"]]

    logger.info(f"Reading unique features from {par['obs_feature']}")
    obs_features = adata.obs[par["obs_feature"]].unique().tolist()

    # format obs_features into file compatible names
    obs_features = [re.sub(r'[\\/*?:"<>|]', "", s) for s in obs_features]
    obs_features = [s.replace(" ", "_") for s in obs_features]
    obs_features = [s.replace("-", "_") for s in obs_features]
    obs_features = [s.replace("'", "") for s in obs_features]

    # generate output dir 
    output_dir = Path(par["output"])
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True)

    logger.info(f"Splitting file based on {par['obs_feature']} values {obs_features}")
    obs_files = []

    for obs_name in obs_features:
        logger.info(f"Filtering modality '{par['modality']}' observations by .obs['{par['obs_feature']}'] == {obs_name}")
        mdata_obs = mdata.copy()

        mdata_obs = mdata_obs[mdata_obs.mod['rna'].obs[par["obs_feature"]] == obs_name]
        mdata_obs_name = f"{input_file.stem}_{obs_name}.h5mu"
        obs_files.append(mdata_obs_name)
        
        logger.info(f"Writing h5mu to file {output_dir / mdata_obs_name}")
        mdata_obs.write_h5mu(output_dir / mdata_obs_name, compression=par["output_compression"])
        
    logger.info(f"Writing output_files CSV file to {par['output_files']}")
    df = pd.DataFrame({"name": obs_features, "filename": obs_files})
    df.to_csv(par["output_files"], index=False)


if __name__ == '__main__':
    main()
