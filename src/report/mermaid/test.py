import unittest
from os import path
import subprocess

class TestMermaid(unittest.TestCase):
    def test_simple_execution(self):
        content = """
            graph LR;
            A--> B & C & D;
            B--> A & E;
            C--> A & E;
            D--> A & E;
            E--> B & C & D;
        """

        input = open("network.input", "w")
        input.write(content)
        input.close()

        out = subprocess.check_output([
              "./mermaid", "--input", "network.input",  "--output", "./network.png"
              ]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("network.png"))

unittest.main()
