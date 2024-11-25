import mudata as mu
import pandas as pd
import re
from pathlib import Path
from collections import defaultdict

### VIASH START
par = {
  'input': 'harmony_knn/integrated.pynndescent_knn.output',
  'modality': 'rna',
  'obs_feature': 'dataset',
  'output': 'reference_download/sample_split',
  'drop_obs_nan': "true",
  'output_compression': None,
  'output_files': 'reference_download/sample_files.csv',
  'ensure_unique_filenames': True
}
import anndata as ad
df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
var3 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
obs3 = pd.DataFrame(["C C", "C_C"], index=df.index, columns=["Obs"])
ad3 = ad.AnnData(df, obs=obs3, var=var3)
mdata = mu.MuData({'rna': ad3})
mdata.write_h5mu("test_san.h5mu")
par["input"] = "test_san.h5mu"
par["obs_feature"] = "Obs"
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

    # sanitize --obs_feature values
    obs_features_s = [re.sub(r'[-\s]', "_", str(s).strip()) for s in obs_features]
    obs_features_s = [re.sub(r'[^A-Za-z0-9_]', "", s) for s in obs_features_s]

    # ensure that names are unique, if not raise or append number as suffix
    if not len(obs_features_s) == len(set(obs_features_s)):
        if not par["ensure_unique_filenames"]:
            raise ValueError(f"File names are not unique after sanitizing the --obs_feature {par['obs_feature']} values")

        logger.info("Ensuring unique names for par['obs_feature']")
        counts = defaultdict(lambda: -1)
        for i, feature in enumerate(obs_features_s):
            counts[feature] += 1
            if (curr_counts := counts[feature]) > 0:
                obs_features_s[i] += f"_{curr_counts}"

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

        # avoid keeping files in memory
        del mdata_obs
        del adata_obs

    logger.info(f"Writing output_files CSV file to {par['output_files']}")
    df = pd.DataFrame({"name": obs_features_s, "filename": obs_files})
    df.to_csv(par["output_files"], index=False)


if __name__ == '__main__':
    main()
