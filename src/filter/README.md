# Filter components

Filter a h5mu file.

Required input arguments:

```yaml
  - name: "--input"
    type: file
    description: Input h5mu file
    direction: input
    required: true
    example: input.h5mu

  - name: "--modality"
    type: string
    multiple: true
    default: "rna"
    required: false
```

Required output arguments:

```yaml
  - name: "--output"
    type: file
    description: Output h5mu file.
    direction: output
    example: output.h5mu

  - name: "--do_subset"
    type: boolean_true
    description: Whether to subset before storing the output.
```

Optional output arguments:

```yaml
  - name: "--obs_name_filter"
    type: string
    default: "filter_<functionality_name>"
    description: In which .obs slot to store a boolean array corresponding to which observations should be removed.

  - name: "--var_name_filter"
    type: string
    default: "filter_<functionality_name>"
    description: In which .var slot to store a boolean array corresponding to which variables should be removed.
```