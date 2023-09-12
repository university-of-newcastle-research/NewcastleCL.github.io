library(rtdists)
library(fmri)
library(pmwg)
library(tibble)
rm(list=ls())

load("fmri.RData")

sampled <- run_stage(sampled, stage = "adapt",iter = 10000, particles = 1000, n_cores = 16, epsilon=0.1)
save.image("fmri.RData")


sampled <- run_stage(sampled, stage = "sample", iter = 1000, particles = 100, n_cores = 16)

save.image("fmri.RData")
