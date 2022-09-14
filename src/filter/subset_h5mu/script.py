import mudata

### VIASH START

par = {
    "input": "foo.h5mu",
    "output": "bar.h5mu",
    "number_of_observations": 2000
}

### VIASH END

if __name__ == "__main__":
    data = mudata.read(par["input"])
    data = data[:par["number_of_observations"]]
    data.write(par["output"], compression="gzip")