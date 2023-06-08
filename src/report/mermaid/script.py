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
# Requests docs on timeout parameter:
# The timeout parameter takes as value a tuple (connect, read)
# The connect timeout is the number of seconds Requests will wait for 
# your client to establish a connection to a remote machine (corresponding to the connect()) 
# call on the socket. It’s a good practice to set connect timeouts to slightly larger than a multiple of 3.
# Once your client has connected to the server and sent the HTTP request, '
# the read timeout is the number of seconds the client will wait for the server to send a response.
response = requests.get('https://mermaid.ink/img/' + base64_string, timeout=(3.05, 10)).content
print(f"Image: {response}", flush=True)
img = Image.open(io.BytesIO(response))
img.save(par["output"],dpi=(600,600))