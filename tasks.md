# Task: Create a new viash component

> Input required: **component name** (e.g. `winnowmap_align`, namespace `winnowmap`)
> Reference: `docs_summary.md` for project conventions and structure.

## Step 1 — Create a git branch

```bash
git checkout -b <component_name>
```

## Step 2 — Create the component directory

```
src/<namespace>/<component_name>/
```

## Step 3 — Create the required files

### 3a. `config.vsh.yaml`
Component metadata and argument definitions. Required fields:
- `name`, `namespace`, `description`, `keywords`
- `links` (homepage, documentation, repository)
- `references` (doi)
- `license`
- `requirements.commands` — list of CLI tools the component calls
- `authors` — use `__merge__` from `/src/_authors/<author>.yaml`
- `argument_groups` — Inputs / Outputs / Options
  - Each argument: `name`, `type`, `description` (markdown), `required`, `example`
  - Outputs need `direction: output`
  - Use `boolean_true` for flags; never add `threads`/`memory` as parameters (use `meta_cpus` / `meta_memory_*` in the script)
- `resources` — `type: bash_script`, `path: script.sh`
- `test_resources` — `test.sh` + `/src/_utils/test_helpers.sh`
- `engines` — Docker image from `quay.io/biocontainers/<tool>:<version>--<build>`; include version detection in setup
- `runners` — `executable` and `nextflow`

### 3b. `script.sh`
Bash wrapper around the tool. Conventions:
- Start with `#!/bin/bash`, then `## VIASH START` / `## VIASH END`, then `set -eo pipefail`
- Unset boolean parameters that are `"false"` before use
- Build arguments as an array `cmd_args=(...)` and call `tool "${cmd_args[@]}"`
- Use `${meta_cpus:+-t "$meta_cpus"}` for threads; never use `par_threads`
- Use 2-space indentation throughout

### 3c. `test.sh`
Unit test script. Conventions:
- Source `"$meta_resources_dir/test_helpers.sh"` and call `setup_test_env`
- Generate test data programmatically with helpers (`create_test_fasta`, `create_test_fastq`, etc.)
- Cover: basic functionality, optional parameters, edge cases
- Validate with `check_file_exists`, `check_file_not_empty`, `check_file_contains`, etc.
- End with `print_test_summary "All tests completed successfully"`

### 3d. `help.txt`
Tool `--help` output for developer reference. No functional purpose.
Format:
```
\```sh
<tool> --help
\```
<paste help output here>
```

### 3e. `test_data/test_data_script.sh`
Script to generate local input files for manual `viash run` testing.
- Use Python3 for reliable data generation (avoids bash pipe issues)
- Generate a reference and query file of realistic size
- Print the `viash run` command at the end

## Step 4 — Update CHANGELOG.md

Add the new component under the current in-progress version block (do **not** create a new version). Format:

```markdown
* `<namespace>`: <one-line description> (PR #xxx):
  - `<namespace>/<component_name>`: <what it does>
```

## Step 5 — Generate test data and run manually

```bash
bash src/<namespace>/<component_name>/test_data/test_data_script.sh

viash run src/<namespace>/<component_name>/config.vsh.yaml -- \
  --<input_arg> src/<namespace>/<component_name>/test_data/<input_file> \
  --output src/<namespace>/<component_name>/test_data/<output_file>
```

Note: always put `--` between `viash run config.vsh.yaml` and the component arguments.

## Step 6 — Run viash tests

```bash
# Test single component
viash test src/<namespace>/<component_name>/config.vsh.yaml

# Test full namespace in parallel
viash ns test --parallel -q <namespace>
```

## Step 7 — Checklist before requesting a review

- [ ] Self-review of all files
- [ ] Conforms to contributing guidelines (see `docs_summary.md`)
- [ ] CHANGELOG.md updated under the current in-progress version (not a new version)
- [ ] `viash ns test --parallel -q <namespace>` passes
- [ ] PR type identified: Breaking change / New functionality / Major / Minor / Docs / Bug fix


---

# Task: Create a Nextflow pipeline using viash components

> Input required: **pipeline name** and the list of **components to chain** (namespace/component_name)
> Reference: `docs_summary.md` for project conventions and structure.

## Background

Each component with `runners: - type: nextflow` in its `config.vsh.yaml` is automatically compiled into a
VDSL3-compatible Nextflow module by `viash ns build`. The generated module lives at:

```
target/nextflow/<namespace>/<component_name>/main.nf
```

It exposes a `run` workflow that follows the standard channel tuple signature:

```
[ id, [ input_map ] ]  →  [ id, [ output_map ] ]
```

## Step 1 — Build all components into NF modules

```bash
viash ns build --setup cb --target target/
```

## Step 2 — Create the pipeline directory

```
src/workflows/<pipeline_name>/
├── main.nf           # Nextflow pipeline script
└── nextflow.config   # executor, container, and parameter defaults
```

No `config.vsh.yaml` is required for a plain Nextflow pipeline.

## Step 3 — Write `main.nf`

```nextflow
nextflow.enable.dsl=2

// Import the required modules (paths relative to pipeline main.nf)
include { <component_name> } from '../../../../target/nextflow/<namespace>/<component_name>/main.nf'

workflow {
  // Build an input channel: [ id, [ input_map ] ]
  Channel.fromPath(params.input)
    | map { f -> [ f.simpleName, [ input: f ] ] }
    | <component_name>.run(
        args: [
          // override component defaults here
        ],
        directives: [
          cpus:   params.cpus   ?: 4,
          memory: params.memory ?: "16 GB"
        ]
      )
    | view { id, outputs -> "Done: $id → ${outputs}" }
}
```

**Chaining multiple components:**

```nextflow
include { component_a } from '.../component_a/main.nf'
include { component_b } from '.../component_b/main.nf'

workflow {
  input_ch
    | component_a.run(args: [...])
    | map { id, out -> [ id, [ input: out.<output_key> ] ] }   // re-map outputs to next inputs
    | component_b.run(args: [...])
}
```

## Step 4 — Write `nextflow.config`

```groovy
docker.enabled = true

params {
  input  = "*.fastq.gz"
  cpus   = 4
  memory = "16 GB"
}

process {
  withName: '<component_name>' {
    container = 'quay.io/biocontainers/<tool>:<version>--<build>'
    cpus      = 4
    memory    = '16 GB'
  }
}
```

The container tag must match exactly what is in the component's `engines.docker.image`.

## Step 5 — Run the pipeline

```bash
nextflow run src/workflows/<pipeline_name>/main.nf \
  --input data/<input_file>
```

## Step 6 — Update CHANGELOG.md

Add the new pipeline under the current in-progress version block:

```markdown
* `workflows`: <one-line description> (PR #xxx):
  - `workflows/<pipeline_name>`: <what it does>
```

## Key conventions

| Concern | Rule |
|---|---|
| Module path | `target/nextflow/<ns>/<comp>/main.nf` — generated, do not edit manually |
| Channel format | Always `[ id, [ map_of_inputs ] ]` — inputs as a named map, never positional |
| Passing args | Use `.run(args: [...])` to override component defaults |
| CPU / memory | Set via `directives` in `.run()` or `process { withName }` in config; never add `par_cpus`/`par_threads` to the component script |
| Re-mapping outputs | Between steps, explicitly map `out.<key>` to the next component's input key |
| Multi-sample | Use `Channel.fromFilePairs` or a samplesheet CSV parsed with `splitCsv` |
