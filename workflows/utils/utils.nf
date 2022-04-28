process emptyCSVOutput {
  tag "${id}"
  echo { (params.debug == true) ? true : false }
  cache 'deep'
  stageInMode "symlink"
  input:
    tuple val(id), path(input), val(pars), val(suffix)
  output:
    tuple val("${id}"), path("${id}.${input}.${suffix}.empty.csv"), val(pars)
  script:
    """
    touch "${id}.${input}.${suffix}.empty.csv"
    """
}


process writeToFile {

  echo { (params.debug == true) ? true : false }
  cache 'deep'
  stageInMode "symlink"
  publishDir "${params.output}", mode: 'copy', overwrite: true
  input:
    tuple val(id), val(input), path(output), val(content)
  output:
    tuple val("${id}"), path("${output}")
  script:
    """
    mkdir -p `dirname $output`
    echo "$content" > $output
    """

}

process rename {

  tag "${id}"
  echo { (params.debug == true) ? true : false }
  cache 'deep'
  stageInMode "symlink"
  input:
    tuple val(id), path(input), val(pars)
  output:
    tuple val("${id}"), path("${id}.${input}"), val(pars)
  script:
    """
    cp $input "${id}.${input}"
    """

}

// A functional approach to 'updating' a value for an option
// in the params Map.
def overrideOptionValue(triplet, _key, _option, _value) {
    mapCopy = triplet[2].toConfigObject().toMap() // As mentioned on https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/config/CascadingConfig.groovy

    return [
        triplet[0],
        triplet[1],
        triplet[2].collectEntries{ function, v1 ->
        (function == _key)
            ? [ (function) : v1.collectEntries{ k2, v2 ->
                (k2 == "arguments")
                    ? [ (k2) : v2.collectEntries{ k3, v3 ->
                        (k3 == _option)
                            ? [ (k3) : v3 + [ "value" : _value ] ]
                            : [ (k3) : v3 ]
                    } ]
                    : [ (k2) : v2 ]
            } ]
            : [ (function), v1 ]
        }
    ]
}


def has_param(name) {
  return params.containsKey(name) && params[name] != ""
}
def check_required_param(name, description) {
  if (!has_param(name)) {
    exit 1, "ERROR: Please provide a --{name} parameter {description}"
  }
}

def getChild(parent, child) {
  if (child.contains("://") || java.nio.file.Paths.get(child).isAbsolute()) {
    child
  } else {
    parent.replaceAll('/[^/]*$', "/") + child
  }
}