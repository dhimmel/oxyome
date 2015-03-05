library(dplyr)
library(ggplot2)

project.dir = '/home/dhimmels/Dropbox/lung/followup/interactome/'

result.df <- file.path(project.dir, 'disease-go-pairs.tsv') %>%
  read.delim()


