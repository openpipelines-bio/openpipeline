library(tidyverse)

lines <- read_lines("https://bitbucket.org/CRSwDev/cwl/raw/8feeace1141b24749ea6003f8e6ad6d3ad5232de/v2.2.1/rhapsody_pipeline_2.2.1.cwl")

js <- jsonlite::fromJSON(lines[-1], simplifyVector = FALSE)

# recursively remove objects with `"class": "DockerRequirement"`
remove_docker <- function(x) {
  if (is.list(x)) {
    if ("class" %in% names(x) && x$class == "DockerRequirement") {
      return(NULL)
    } else {
      out <- map(x, remove_docker)
      out2 <- out[!sapply(out, is.null)]
      return(out2)
    }
  } else {
    return(x)
  }
}
js2 <- remove_docker(js)
# js2 <- js

lines2 <- c(lines[[1]], jsonlite::toJSON(js2, auto_unbox = TRUE, pretty = TRUE))

write_lines(lines2, "src/mapping/bd_rhapsody_2/rhapsody_pipeline_2.2.1_nodocker.cwl")
