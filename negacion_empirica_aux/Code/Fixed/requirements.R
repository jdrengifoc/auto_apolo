packages <- c(
  "doParallel", "snow", "fdrtool", "plm"
)

if (!all(packages %in% rownames(installed.packages()))){
  idx <- which((packages %in% rownames(installed.packages()))==F)
  install.packages(packages[idx])  
  rm(idx)
}

lapply(packages, library, character.only = TRUE) 
rm(packages)
