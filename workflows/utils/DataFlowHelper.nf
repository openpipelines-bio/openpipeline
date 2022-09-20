/* usage:
| splitParams(
  pca: [ "input": "input", "obsm_output": "obsm_pca" ]
  harmonypy: [ "obs_covariates": "obs_covariates", "obsm_input": "obsm_pca" ],
  find_neighbors: [ "obsm_input": "obsm_pca" ],
  umap: [ "output": "output" ]
)
*/

def splitParams(Map args) {
  wfKey = args.key ?: "splitParams"
  args.keySet().removeAll(["key"])

  
  /*
  data = [a:1, b:2, c:3]
  // args = [foo: ["a", "b"], bar: ["b"]]
  args = [foo: [a: 'a', out: "b"], bar: [in: "b"]]
  */
  
  workflow splitParamsInstance {
    take:
    input_

    main:
    output_ = input_
      | map{ tup -> 
        id = tup[0]
        data = tup[1]
        passthrough = tup.drop(2)

        // determine new data
        toRemove = args.collectMany{ _, dataKeys -> 
          // dataKeys is a map but could also be a list
          dataKeys instanceof List ? dataKeys : dataKeys.values()
        }.unique()
        newData = data.findAll{!toRemove.contains(it.key)}

        // determine splitargs
        splitArgs = args.collectEntries{procKey, dataKeys -> 
          // dataKeys is a map but could also be a list
          newSplitData = dataKeys.collectEntries{ val ->
            newKey = val instanceof String ? val : val.key
            origKey = val instanceof String ? val : val.value
            [ newKey, data[origKey] ]
          }
          [procKey, newSplitData]
        }

        // return output
        [ id, newData, splitArgs] + passthrough
      }

    emit:
    output_
  }

  return splitParamsInstance.cloneWithName(wfKey)
}

/* usage:
| combineParams("harmonypy")
*/


def combineParams(String key) {
  wfKey = "combineParams_" + key
  
  workflow combineParamsInstance {
    take:
    input_

    main:
    output_ = input_
      | map{ tup -> 
        id = tup[0]
        data = tup[1]
        splitArgs = tup[2]
        passthrough = tup.drop(3)

        // try to infer arg name
        if (data !instanceof Map) {
          data = [ input: data ]
        }
        newData = data + splitArgs[key]

        [ id, newData, splitArgs] + passthrough
      }

    emit:
    output_
  }

  return combineParamsInstance.cloneWithName(wfKey)

}

