library(tidyverse)

workflows <- yaml::yaml.load(system("viash ns list -q '^workflows'", intern = TRUE))

outs <- map_df(workflows, function(wf) {
  cat("Running ", wf$functionality$namespace, "/", wf$functionality$name, "\n", sep = "")

  dir <- dirname(wf$info$config)
  tests <- wf$functionality$test_resources %>% keep(~ .$type == "nextflow_script")

  if (length(tests) == 0) {
    tibble(
      namespace = wf$functionality$namespace,
      functionality = wf$functionality$name,
      platform = "native",
      test_name = "tests",
      exit_code = -1L,
      duration = 0L,
      result = "MISSING"
    )
  } else {
    map_df(
      tests,
      function(test) {
        if (file.exists(paste0(dir, "/graph.dot"))) file.remove(paste0(dir, "/graph.dot"))
        
        args <- c(
          "run", ".",
          "-main-script", paste0(dir, "/", test$path),
          "-entry", test$entrypoint,
          "-profile", "docker,no_publish",
          "-with-dag", paste0(dir, "/graph.dot"),
          "-resume"
        )

        start_time <- Sys.time()
        out <- processx::run(
          "nextflow",
          args = args,
          error_on_status = FALSE,
          env = c("current", NXF_VER = "22.04.5")
        )
        stop_time <- Sys.time()
        duration <- ceiling(as.numeric(difftime(stop_time, start_time, unit = "sec")))
        result <- if (out$status > 0) "FAILED" else "SUCCESS"
        if (out$status > 0) {
          cat(
            "========================== ERROR LOG ==========================\n",
            out$stdout,
            "===============================================================\n",
            sep = ""
          )
        }

        tibble(
          namespace = wf$functionality$namespace,
          functionality = wf$functionality$name,
          platform = "native",
          test_name = paste0(basename(test$path), "$", test$entrypoint),
          exit_code = out$status,
          duration = duration,
          result = result
        )
      }
    )
  }
})

write_tsv(outs, ".viash_log_integration.tsv")