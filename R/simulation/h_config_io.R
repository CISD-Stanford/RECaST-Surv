# Explicit data-only configuration loading, separate from the pure DGM.
# package_root is this isolated package, never the main campaign or its runtime.
H_read_configs <- function(package_root, scenario) {
  H_require(is.character(scenario) && length(scenario) == 1L &&
    scenario %in% c("M1", paste0("H", 1:5)), "Unknown configuration document")
  configs <- dget(file.path(package_root, "config", paste0(scenario, ".R")))
  H_require(is.list(configs) && length(configs) > 0L, "Empty configuration document")
  for (cfg in configs) {
    H_validate_config(cfg)
    H_require(identical(cfg$scenario, scenario), "Wrong scenario document")
  }
  variants <- vapply(configs, function(cfg) cfg$variant, character(1))
  expected <- if (scenario == "H1") as.vector(outer(c(.5, 1, 2), c(0, -.5),
    function(m, d) paste0("m", m, "_Delta", d))) else if (scenario == "H3")
      c("all_zero_null", "delayed_alternative") else paste0("Delta", c(0, -.5))
  H_require(identical(sort(variants), sort(expected)), "Incomplete, duplicate or unapproved variants")
  configs
}
