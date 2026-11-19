"""Run cwltool, limiting the resources used by its parallel executor.

cwltool's parallel executor decides which steps can run at the same time based
on the memory and cpus of the host. When running in a container, the resources
allocated to the container are usually lower than those of the host.
"""

import argparse
import sys

from cwltool.executors import MultithreadedJobExecutor
from cwltool.main import main

parser = argparse.ArgumentParser(allow_abbrev=False)
parser.add_argument("--max_ram_mb", type=int)
parser.add_argument("--max_cores", type=int)
args, cwltool_args = parser.parse_known_args()

executor = None
if "--parallel" in cwltool_args:
    executor = MultithreadedJobExecutor()
    # These attributes are set from the resources of the host when creating the executor.
    # Fail when they do not exist, instead of silently using the resources of the host.
    for attribute, value in [
        ("max_ram", args.max_ram_mb),
        ("max_cores", args.max_cores),
    ]:
        if not hasattr(executor, attribute):
            raise AttributeError(
                f"cwltool's MultithreadedJobExecutor has no attribute '{attribute}'"
            )
        if value:
            setattr(executor, attribute, value)

sys.exit(main(argsl=cwltool_args, executor=executor))
