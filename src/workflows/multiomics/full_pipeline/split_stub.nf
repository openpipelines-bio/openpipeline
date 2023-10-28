process splitStub {
      input:
        tuple val(id), path(unused), val(passthrough)

      output:
        tuple val(id), path("stub_h5mus"), path("modalities.csv"), val(passthrough)

      script:
        """
        echo "This process is not meant to be run without -stub being defined."
        exit 1
        """

      stub:
        """
        mkdir stub_h5mus
        touch stub_h5mus/${id}_rna.h5mu
        touch stub_h5mus/${id}_prot.h5mu
        touch stub_h5mus/${id}_vdj.h5mu
        echo -e "name,filename\nrna,stub_h5mus/${id}_rna.h5mu\nprot,stub_h5mus/${id}_prot.h5mu\nvdj,stub_h5mus/${id}_vdj.h5mu" > modalities.csv
        """
    }