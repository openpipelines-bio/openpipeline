def pytest_configure(config):
    config.addinivalue_line(
        "markers", "modality: only run tests that target specific modality"
    )
    config.addinivalue_line(
        "markers", "layer: only run tests that target specified layer"
    )
