
library(tidyverse)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)

infile <- read.csv(args[1],sep='\t') 
out <- infile %>%
    select(-c(Sample,Chromosome, Stop_Position, Potential_sequence, Alignment)) %>%
    filter(Potential_IS != "No identity") %>%
    pivot_wider(names_from = Potential_IS, values_from = Start_Position, values_fn = list)
out <- apply(out,2,as.character)

write.table(out, args[2], row.names=FALSE, sep="\t",quote=FALSE)