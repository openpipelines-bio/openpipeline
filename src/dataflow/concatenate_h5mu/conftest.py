import inspect

import pytest
from pytest_benchmark.fixture import BenchmarkFixture

# Workaround: pytest-benchmark reads the private 'node._nodeid', which was removed in pytest 9.
# Fails the test run as soon as pytest-benchmark is fixed, so that this file can be removed.


def pytest_collection_modifyitems(items):
    benchmark_source = inspect.getsource(BenchmarkFixture.__init__)
    if "node.nodeid" in benchmark_source or "node._nodeid" not in benchmark_source:
        raise pytest.UsageError(
            f"pytest-benchmark no longer reads 'node._nodeid', "
            f"the workaround in {__file__} is no longer needed. Please remove it."
        )
    for item in items:
        item._nodeid = item.nodeid
