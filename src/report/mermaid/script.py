import base64
import requests, io
from PIL import Image

## VIASH START
par = {
}
## VIASH END

file = open(par["input"],mode='r')
graph = file.read()
file.close()

graphbytes = graph.encode("ascii")
base64_bytes = base64.b64encode(graphbytes)
base64_string = base64_bytes.decode("ascii")
print(f"Ascii: {base64_string}", flush=True)
response = io.BytesIO(requests.get('https://mermaid.ink/img/' + base64_string).content)
print(f"Image: {response}", flush=True)
img = Image.open(response)
img.save(par["output"],dpi=(600,600))