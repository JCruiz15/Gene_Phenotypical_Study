pkg_list <- c("BiocManager", "data.table")
bioc_list <- c("STRINGdb", "igraph", "linkcomm")

for (pkg in pkg_list) {
  install.packages(pkg, repos = "https://cran.rediris.es/", 
                   lib = "./R_deps") 
  if (!require(pkg, character.only = T)) {
    quit(save = "no", status = 1, runLast = FALSE)
  }
}

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

for (biopkg in bioc_list) {
  BiocManager::install(biopkg, lib = "./R_deps", force=TRUE)
}
