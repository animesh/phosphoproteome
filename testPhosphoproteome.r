source("R/Phosphoproteome.R")
library()
install.packages("magrittr") # package installations are only needed the first time you use %>% as per https://stackoverflow.com/questions/30248583/error-could-not-find-function
library(magrittr)
install.packages("stringr") # for str_split
library(stringr)
data<-PhosphorylatedProteinInformation("201105_kath_phostot_1B-(prots).xlsx")
