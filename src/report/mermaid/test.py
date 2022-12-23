import pytest
from os import path
import sys

content = """\
graph LR;
    A--> B & C & D;
    B--> A & E;
    C--> A & E;
    D--> A & E;
    E--> B & C & D;
"""

def test_simple_execution(run_component):

    with open("network.input", "w") as file:
        file.write(content)

    run_component([
        "--input", "network.input",
        "--output", "network.png"
    ])
    
    assert path.exists("network.png")

if __name__ == '__main__':
    sys.exit(pytest.main([__file__], plugins=["viashpy"]))
