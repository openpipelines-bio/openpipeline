library(tidyverse)

# url <- "https://raw.githubusercontent.com/alexdobin/STAR/master/source/parametersDefault"
# param_txt <- readr::read_lines(url)

# until https://github.com/alexdobin/STAR/pull/1710 is merged
param_txt <- readr::read_lines("src/mapping/star_align/utils/parametersDefault")

dev_begin <- grep("#####UnderDevelopment_begin", param_txt)
dev_end <- grep("#####UnderDevelopment_end", param_txt)

# strip development sections
nondev_ix <- unlist(map2(c(1, dev_end + 1), c(dev_begin - 1, length(param_txt)), function(i, j) {
  if (i >= 1 && i < j) {
    seq(i, j, 1)
  } else {
    NULL
  }
}))

param_txt2 <- param_txt[nondev_ix]

# strip comments
param_txt3 <- param_txt2[-grep("^#[^#]", param_txt2)]

# detect groups
group_ix <- grep("^### ", param_txt3)

out <- map2_dfr(
  group_ix,
  c(group_ix[-1] - 1, length(param_txt3)),
  function(group_start, group_end) {
    # cat("group_start <- ", group_start, "; group_end <- ", group_end, "\n", sep = "")
    group_name <- gsub("^### ", "", param_txt3[[group_start]])

    group_txt <- param_txt3[seq(group_start + 1, group_end)]

    arg_ix <- grep("^[^ ]", group_txt)

    arguments <- map2_dfr(
      arg_ix,
      c(arg_ix[-1] - 1, length(group_txt)),
      function(arg_start, arg_end) {
        # cat("arg_start <- ", arg_start, "; arg_end <- ", arg_end, "\n", sep = "")

        # process name and default
        first_txt <- group_txt[[arg_start]]
        first_regex <- "^([^ ]*) +(.*) *$"
        if (!grepl(first_regex, first_txt)) {
          stop("Line '", first_txt, "' did not match regex '", first_regex, "'")
        }
        name <- gsub(first_regex, "\\1", first_txt)
        default <- gsub(first_regex, "\\2", first_txt)

        # process type and first description
        second_txt <- group_txt[[arg_start + 1]]
        second_regex <- "^   +([^:]*):[ ]+(.*)$"
        if (!grepl(second_regex, second_txt)) {
          stop("Line '", second_txt, "' did not match regex '", second_regex, "'")
        }
        type <- gsub(second_regex, "\\1", second_txt)
        desc_start <- str_trim(gsub(second_regex, "\\2", second_txt))

        # process more description
        desc_cont1 <- group_txt[seq(arg_start + 2, arg_end)]

        desc <-
          if (sum(str_length(desc_cont1)) == 0) {
            desc_start
          } else {
            # detect margin
            margins <- str_extract(desc_cont1, "^( +)") %>% na.omit
            margin <- margins[[which.min(str_length(margins))]]
            desc_cont2 <- gsub(paste0("^", margin), "", desc_cont1)
            desc_cont3 <- ifelse(grepl("\\.\\.\\.", desc_cont2), paste0("- ", desc_cont2), desc_cont2)
            desc_cont4 <- str_trim(desc_cont3)

            # construct desc
            str_trim(paste0(c(desc_start, "", desc_cont4), "\n", collapse = ""))
          }

        tibble(
          group_name,
          name,
          default,
          type,
          description = desc
        )
      }
    )

    arguments
  }
)

# todo: manually fix alignEndsProtrude?
# assigning types
type_map <- c("string" = "string", "int" = "integer", "real" = "double", "double" = "double", "int, string" = "string")
file_args <- c("genomeDir", "readFilesIn", "sjdbGTFfile", "genomeFastaFiles", "genomeChainFiles", "readFilesManifest")
long_args <- c("limitGenomeGenerateRAM", "limitIObufferSize", "limitOutSAMoneReadBytes", "limitBAMsortRAM")
required_args <- c("genomeDir", "readFilesIn")

# converting examples
as_safe_int <- function(x) tryCatch({as.integer(x)}, warning = function(e) { bit64::as.integer64(x) })
safe_split <- function(x) strsplit(x, "'[^']*'(*SKIP)(*F)|\"[^\"]*\"(*SKIP)(*F)|\\s+", perl = TRUE)[[1]] %>% gsub("^[\"']|[\"']$", "", .)
trafos <- list(
  string = function(x) x,
  integer = as_safe_int,
  double = as.numeric,
  strings = function(x) safe_split(x),
  integers = function(x) sapply(safe_split(x), as_safe_int),
  doubles = function(x) as.numeric(safe_split(x))
)
# remove arguments that are not relevant for viash
removed_args <- c("versionGenome", "parametersFiles", "sysShell", "runDirPerm")
# these settings are defined by the viash component
manual_args <- c("runThreadN", "outTmpDir", "runMode", "outFileNamePrefix", "genomeDir", "readFilesIn")

# make viash-like values
out2 <- out %>%
  # remove arguments that are not relevant for viash
  filter(!name %in% c(removed_args, manual_args)) %>%
  # remove arguments that are related to a different runmode
  filter(!grepl("--runMode", description) | grepl("--runMode alignReads", description)) %>%
  filter(!grepl("--runMode", group_name) | grepl("--runMode alignReads", group_name)) %>%
  mutate(
    viash_arg = paste0("--", name),
    type_step1 = type %>%
      str_replace_all(".*(int, string|string|int|real|double)\\(?(s?).*", "\\1\\2"),
    viash_type = type_map[gsub("(int, string|string|int|real|double).*", "\\1", type_step1)],
    multiple = type_step1 == "int, string" | grepl("s$", type_step1) | grepl("^[4N][\\* ]", type),
    default_step1 = default %>%
      {ifelse(. %in% c("-", "None"), NA_character_, .)},
    viash_default =
      mapply(
        default_step1,
        paste0(viash_type, ifelse(multiple, "s", "")),
        FUN = function(str, typ) trafos[[typ]](str)
      ),
    # viash_type = ifelse(sapply(viash_default, bit64::is.integer64), "long", viash_type),
    # update type
    viash_type = case_when(
      name %in% long_args ~ "long",
      name %in% file_args ~ "file",
      TRUE ~ viash_type
    ),
    # turn longs into character because yaml::write_yaml doesn't handle longs well
    viash_default = ifelse(sapply(viash_default, bit64::is.integer64), map(viash_default, as.character), viash_default),
    group_name = gsub(" - .*", "", group_name),
    required = ifelse(name %in% required_args, TRUE, NA)
  )
print(out2, n = 200)
out2 %>% mutate(i = row_number()) %>% 
  # filter(is.na(default_step1) != is.na(viash_default)) %>%
  select(-group_name, -description)

out2 %>% filter(!grepl("--runMode", description) | grepl("--runMode alignReads", description))

argument_groups <- map(unique(out2$group_name), function(group_name) {
  args <- out2 %>%
    filter(group_name == !!group_name) %>%
    pmap(function(viash_arg, viash_type, multiple, viash_default, description, required, ...) {
      li <- lst(
        name = viash_arg,
        type = viash_type,
        description = description
      )
      if (all(!is.na(viash_default))) {
        li$example <- viash_default
      }
      if (!is.na(multiple) && multiple) {
        li$multiple <- multiple
        li$multiple_sep <- ";"
      }
      if (!is.na(required) && required) {
        li$required <- required
      }
      li
    })
  list(name = group_name, arguments = args)
})

yaml::write_yaml(list(argument_groups = argument_groups), "src/mapping/star_align/argument_groups.yaml")
