library(tidyverse)

# fetch and read cwl file
lines <- read_lines("https://bitbucket.org/CRSwDev/cwl/raw/8feeace1141b24749ea6003f8e6ad6d3ad5232de/v2.2.1/rhapsody_pipeline_2.2.1.cwl")
cwl_header <- lines[[1]]
cwl_obj <- jsonlite::fromJSON(lines[-1], simplifyVector = FALSE)

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
cwl_obj_2 <- remove_docker(cwl_obj)
# js2 <- js

# write to file
new_lines <- c(cwl_header, jsonlite::toJSON(cwl_obj_2, auto_unbox = TRUE, pretty = TRUE))

write_lines(new_lines, "src/bd_rhapsody/bd_rhapsody_sequence_analysis/rhapsody_pipeline_2.2.1_nodocker.cwl")

# detect main workflow arguments
gr <- dynutils::list_as_tibble(cwl_obj_2$`$graph`)

gr %>% print(n = 100)

main <- gr %>% filter(gr$id == "#main")

main_inputs <- main$inputs[[1]]

input_ids <- main_inputs %>% map_chr("id") %>% gsub("^#main/", "", .)

# check whether in config
config <- yaml::read_yaml("src/bd_rhapsody/bd_rhapsody_sequence_analysis/config.vsh.yaml")
config$all_arguments <- config$argument_groups %>% map("arguments") %>% list_flatten()
arg_names <- config$all_arguments %>% map_chr("name") %>% gsub("^--", "", .)

# arguments in cwl but not in config
setdiff(tolower(input_ids), arg_names)

# arguments in config but not in cwl
setdiff(arg_names, tolower(input_ids))

# create arguments from main_inputs
arguments <- map(main_inputs, function(main_input) {
  input_id <- main_input$id %>% gsub("^#main/", "", .)
  input_type <- main_input$type[[2]]

  if (is.list(input_type) && input_type$type == "array") {
    multiple <- TRUE
    input_type <- input_type$items
  } else {
    multiple <- FALSE
  }

  if (is.list(input_type) && input_type$type == "enum") {
    choices <- input_type$symbols %>%
      gsub(paste0(input_type$name, "/"), "", .)
    input_type <- "enum"
  } else {
    choices <- NULL
  }

  description <-
    if (is.null(main_input$label)) {
      main_input$doc
    } else if (is.null(main_input$doc)) {
      main_input$label
    } else {
      paste0(main_input$label, ". ", main_input$doc)
    }

  type_map <- c(
    "float" = "double",
    "int" = "integer",
    "string" = "string",
    "boolean" = "boolean",
    "File" = "file",
    "enum" = "string"
  )

  out <- list(
    name = paste0("--", tolower(input_id)),
    type = type_map[input_type],
    # TODO: use summary when viash 0.9 is released
    # summary = main_input$doc,
    # description = main_input$doc,
    description = description,
    multiple = multiple,
    choices = choices,
    info = list(
      config_key = input_id
    )
  )

  out[!sapply(out, is.null)]
})



yaml::write_yaml(
  arguments,
  "src/bd_rhapsody/bd_rhapsody_sequence_analysis/arguments.yaml",
  handlers = list(
    logical = yaml::verbatim_logical
  )
)
