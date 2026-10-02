import json
import os
import shutil
import subprocess
import sys
import tempfile

from plotly.offline import get_plotlyjs

## VIASH START
par = {
    "input": "src/report/deseq2_report/test_data/overall",
    "input_prefix": "deseq2_analysis",
    "output": "report.html",
    "output_data": None,
    "title": "Differential expression report",
    "project": None,
    "design_description": None,
    "pipeline": None,
    "reference": None,
    "attribution": None,
    "counts_description": None,
    "obs_sample_label": None,
    "group_labels": None,
    "obs_pair": None,
    "pair_label": None,
    "highlight_genes": None,
    "p_adj_threshold": None,
    "log2fc_threshold": None,
    "min_count": 10,
    "n_pca_genes": 500,
    "n_heatmap_genes": 28,
    "n_table_genes": 300,
    "max_volcano_ns_genes": 2500,
    "n_highlight_genes": 4,
    "min_basemean_highlight": 100,
    "pvalue_bins": 40,
    "seed": 0,
}
meta = {"resources_dir": "src/report/deseq2_report", "temp_dir": "/tmp"}
## VIASH END

sys.path.append(meta["resources_dir"])
from report_data import build_report_data  # noqa: E402
from setup_logger import setup_logger  # noqa: E402

logger = setup_logger()

# Files rendered together with report.qmd
REPORT_FILES = [
    "report.qmd",
    "theme-light.scss",
    "theme-dark.scss",
    "report.css",
    "logo.html",
]


def render_report(data, output):
    """Render report.qmd with the report data into a self-contained HTML file."""
    with tempfile.TemporaryDirectory(dir=meta["temp_dir"]) as tmp:
        for name in REPORT_FILES:
            shutil.copy(os.path.join(meta["resources_dir"], name), tmp)
        # plotly.js once in the header, so the report works offline
        with open(os.path.join(tmp, "plotly.html"), "w") as f:
            f.write(f"<script>{get_plotlyjs()}</script>\n")
        # report.qmd reads the data from report_data.json next to it
        with open(os.path.join(tmp, "report_data.json"), "w") as f:
            json.dump(data, f, allow_nan=False)

        # Quarto, Deno and Jupyter write caches and runtime files to the home
        # directory, which is not always writable in a container
        env = dict(
            os.environ,
            HOME=tmp,
            XDG_CACHE_HOME=os.path.join(tmp, ".cache"),
            XDG_DATA_HOME=os.path.join(tmp, ".local", "share"),
            DENO_DIR=os.path.join(tmp, ".deno"),
            JUPYTER_RUNTIME_DIR=os.path.join(tmp, ".jupyter"),
            IPYTHONDIR=os.path.join(tmp, ".ipython"),
            QUARTO_PYTHON=sys.executable,
        )
        meta_text = data["meta"]
        subtitle_text = "  -  ".join(
            x for x in [meta_text["project"], meta_text["design"]] if x
        )
        subtitle = ["-M", f"subtitle:{subtitle_text}"] if subtitle_text else []
        logger.info("Rendering the report with Quarto")
        subprocess.run(
            [
                "quarto",
                "render",
                "report.qmd",
                "--output",
                "report.html",
                "-M",
                f"title:{data['meta']['title']}",
            ]
            + subtitle,
            cwd=tmp,
            env=env,
            check=True,
        )
        shutil.copy(os.path.join(tmp, "report.html"), output)


def main():
    data = build_report_data(par, logger)
    if par["output_data"]:
        with open(par["output_data"], "w") as f:
            json.dump(data, f, allow_nan=False)
        logger.info("Wrote report data to %s", par["output_data"])
    render_report(data, par["output"])
    logger.info("Wrote report to %s", par["output"])


if __name__ == "__main__":
    main()
