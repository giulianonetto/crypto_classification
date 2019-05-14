suppressPackageStartupMessages({
  suppressWarnings({
    library(tidyverse)
    library(ggpubr)
    library(reshape2)
    library(EBImage)
    library(heatmaply)
  })
})

# set global ggplot theme
theme_set(theme_pubclean(base_size = 16))

# load data
testset = read.csv("modelImages/preprocessed/testing_set.tsv", 
                   sep = "\t", header = TRUE)
trainset = read.csv("modelImages/preprocessed/training_set.tsv", 
                    sep = "\t", header = TRUE)
big.list = readRDS("big.list.RDS")
fit = readRDS("tuning_fit_random_forest.RDS") # it has all caret's objects
pred_probs = readRDS("predicted_probabilities.RDS")
preds = readRDS("predicted_classes.RDS")
df = data.frame("Predicted" = preds, "Reference" = testset$Type,
                pred_probs)
rownames(df) = NULL
filepath = "modelImages/plots/fig6_heatmap_probabilities.jpeg"
h = heatmaply(df, file = filepath, 
              showticklabels = c(T, F), 
              width = 700, height = 750, fontsize_col = 14)
htmlwidgets::saveWidget(h,"heatmap_probabilities.html")
system("mv heatmap_probabilities.html modelImages/plots")
