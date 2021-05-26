nextflow.enable.dsl=2

println ""
println "# debug info 'workflow'"
println "$workflow"

println ""
println "# debug info 'nextflow'"
println "$nextflow"

println ""

println "${workflow.projectDir}"
println "${projectDir}"
