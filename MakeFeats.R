source("DataFeatures.R")

data = read.csv("testdata.csv") #example dataset

feats <- Allfeatures(data) #Takes a dataset normalised between 0 and 1

write.csv(feats,"feats.csv")