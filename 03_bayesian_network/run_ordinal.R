if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
project_root <- getwd()
for (file in c("00_packages.R", "01_config.R", "02_data.R", "03_model.R", "05_networks.R", "06_metrics.R", "08_tables.R", "09_ordinal.R"))
  source(file, encoding = "UTF-8")
args <- commandArgs(TRUE)
action <- if (length(args)) args[1] else "plan"
selector <- if (length(args) > 1L) args[2] else "all"
if (!action %in% c("plan", "fit", "summarize")) stop("Usage: run_ordinal.R plan|fit|summarize [all|support|anxiety|mpiq|node1,node2]")
ctx <- ordinal_context()
if (action == "summarize") summarize_ordinal(ctx) else fit_ordinal_targets(ctx, selector, action)
