../bin/viash build ../src/report/mermaid/config.vsh.yaml -o ../target/report/mermaid/

../target/report/mermaid/mermaid -i pipelines-target-p1.mermaid -o pipelines-target-p1.png
../target/report/mermaid/mermaid -i pipelines-target-p2.mermaid -o pipelines-target-p2.png
../target/report/mermaid/mermaid -i pipelines-target-p3.mermaid -o pipelines-target-p3.png