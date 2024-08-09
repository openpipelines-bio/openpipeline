import mudata as mu
import pandas as pd
import re
from pathlib import Path

### VIASH START
par = {
  'input': 'harmony_knn/integrated.pynndescent_knn.output',
  'modality': 'rna',
  'obs_feature': 'dataset',
  'output': 'reference_download/sample_split',
  'drop_obs_nan': "true",
  'output_compression': None,
  'output_files': 'reference_download/sample_files.csv'
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
    obs_features_s = [s.strip().replace(" ", "_") for s in obs_features]
    obs_features_s = [re.sub(r'[-\s]', "_", s) for s in obs_features]
    obs_features_s = [re.sub(r'[^A-Za-z0-9_]', "", s) for s in obs_features]

    # ensure that names are unique, if not append number as suffix
    if not len(obs_features_s) == len(set(obs_features_s)):
        logger.info("Ensuring unique names for par['obs_feature']")
        counts = {}
        for i, feature in enumerate(obs_features_s):
            if feature not in counts:
                counts[feature] = 0
            else:
                counts[feature] += 1
                obs_features_s[i] += "_" + str(counts[feature])

    # generate output dir
    output_dir = Path(par["output"])
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True)

    # split modality of mdata file base on obs_feature
    logger.info(f"Splitting file based on {par['obs_feature']} values {obs_features}")
    obs_files = []

    for obs_name, file_name in zip(obs_features, obs_features_s):
        logger.info(f"Filtering modality '{par['modality']}' observations by .obs['{par['obs_feature']}'] == {obs_name}")
        mdata_obs = mdata.copy()
        adata_obs = mdata_obs.mod[par["modality"]]

        # split the samples
        adata_obs = adata_obs[adata_obs.obs[par["obs_feature"]] == obs_name]
        mdata_obs_name = f"{input_file.stem}_{file_name}.h5mu"
        obs_files.append(mdata_obs_name)

        # Dropping columns that only have nan values after splitting
        if par["drop_obs_nan"]:
            logger.info(f"Dropping all .obs columns with NaN values")
            adata_obs.obs.dropna(axis=1, how='all', inplace=True)

        # replace mdata file with modality adata contianing split samples
        logger.info(f"Writing h5mu filtered for {par['obs_feature']} {obs_name} to file {output_dir / mdata_obs_name}")
        mdata_obs.mod[par["modality"]] = adata_obs
        mdata_obs.write_h5mu(output_dir / mdata_obs_name, compression=par["output_compression"])

    logger.info(f"Writing output_files CSV file to {par['output_files']}")
    df = pd.DataFrame({"name": obs_features_s, "filename": obs_files})
    df.to_csv(par["output_files"], index=False)


if __name__ == '__main__':
    main()
