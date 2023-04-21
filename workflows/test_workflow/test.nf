process testProcess {
      input:
        tuple val(id), path(input_file)
    
      output:
        tuple val(id), path("test_output")

      script:
        """
        cat ${input_file} > test_output
        """
    }