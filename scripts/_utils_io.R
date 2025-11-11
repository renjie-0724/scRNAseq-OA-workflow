suppressPackageStartupMessages({
  library(optparse); library(yaml); library(stringr); library(glue)
})

read_cfg <- function(cfg){
  y <- yaml::read_yaml(cfg)
  # 
  dir.create("logs", showWarnings = FALSE, recursive = TRUE)
  dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
  dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
  y
}

opt_config <- function(){
  OptionParser(option_list = list(
    make_option(c("-c","--config"), type="character", default="config/config.yaml",
                help="Path to config YAML")
  )) |> parse_args() |> (\(o) o$config)()
}

saveRDS_safe <- function(obj, path){
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
}
