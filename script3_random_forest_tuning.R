suppressPackageStartupMessages({
  suppressWarnings({
    library(caret)
    library(caTools)
    library(ggpubr)
    library(tidyverse)
  })
})


# set theme globally (ggpubr for publication-ready plots)
theme_set(theme_pubr())

# Dataset split
print(str_glue("Loading and preprocessing initial data"))

df = read.csv("features_all_labeled.tsv", 
              sep = "\t", header = TRUE, stringsAsFactors = FALSE)
set.seed(829)
split <- sample.split(df$Type, SplitRatio = 0.75)
trainset = df[split,]
testset = df[!split,]
test_cell_ids = testset$Cell
train_cell_ids = trainset$Cell
idslist = list("test_ids" = test_cell_ids, "train_ids" = train_cell_ids)
saveRDS(idslist, "cell_ids.RDS")
trainset$Cell = NULL
testset$Cell = NULL
normParams = preProcess(trainset, 
                        method = c("BoxCox", "center", "scale"))
testset = predict(normParams, testset)
trainset = predict(normParams, trainset)

# save training and testing sets
print(str_glue("Saving training and testing sets files"))
write.table(trainset,
            "modelImages/preprocessed/training_set.tsv",
            sep = "\t", quote = F, row.names = T,
            col.names = T)
write.table(testset,
            "modelImages/preprocessed/testing_set.tsv",
            sep = "\t", quote = F, row.names = T,
            col.names = T)

# Model tuning
rf_ctrl = trainControl(method = "repeatedcv",
		       number = 10,
                       ## repeated ten times
                       repeats = 5,
		       classProbs = TRUE,
                       savePredictions = TRUE,
                       summaryFunction = multiClassSummary)
rf_grid = expand.grid(minNode = c(1:10),
                      predFixed = c(4:20))

print(str_glue("Tuning Random Forest classifier parameters - this may take awhile..."))

# WARNING: takes forever!!
if (!file.exists("tuning_fit_random_forest.RDS")) {
  set.seed(5996)
  fit_rf = caret::train(Type ~ ., data = trainset,
                        method = "Rborist",
                        tuneGrid = rf_grid,
                        trControl = rf_ctrl,
                        nTree = 2000,
                        nthread = 4)
  saveRDS(fit_rf, "tuning_fit_random_forest.RDS")
} else {
  fit_rf = readRDS("tuning_fit_random_forest.RDS")
}

print(str_glue("Done with model tuning!"))

best_predFixed = fit_rf$best$predFixed
best_minNode = fit_rf$bestTune$minNode
accuracy = max(fit_rf$results$Accuracy) * 100

print(str_glue("Best number of randomly selected predictors: {best_predFixed}"))
print(str_glue("Best minimum node size: {best_minNode}"))
print(str_glue("Classifier's top overall accuracy: {round(accuracy, 2)} %"))

# build tuning plot
tuning_plot_rf = ggplot(fit_rf, highlight = TRUE) +
  ggtitle("Random Forest - Training Accuracy (nTree = 500)") +
  theme(plot.title = element_text(hjust = 0.475, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(hjust = 0, face = "bold")) +
  labs(color = "Minimum\nNode Size", shape = "Minimum\nNode Size") +
  scale_shape_manual(values = 1:10) + 
  scale_color_manual(values = 1:10) +
  labs(y = "Accuracy (Repeated CV)")

# update plot's scale and annotations
tuning_plot_rf$data$Accuracy = tuning_plot_rf$data$Accuracy * 100
big = max(tuning_plot_rf$data$Accuracy) %>% round(2)
small = min(tuning_plot_rf$data$Accuracy) 
plot_text = str_glue("{best_minNode} Nodes Minimum
                     {best_predFixed} Randomly Selected Predictors
                     Highest Training Accuracy = {big} %")
tuning_plot_rf$labels$y = str_glue("{tuning_plot_rf$labels$y} %")
tuning_plot_rf = tuning_plot_rf + 
  ylim(small, big*1.006) +
  geom_hline(yintercept = big, linetype = "longdash", color = "gray40") +
  annotate("text", x = 10, y = big * 1.0025, hjust = 0,
           fontface = "bold", color = "gray10",
           label = plot_text)

# save plot
plot_path = "modelImages/plots/fig2_tuning_random_forest.jpeg"
print(str_glue("Saving tuning plot in:\n {plot_path}"))
ggsave(plot_path, tuning_plot_rf, width = 7.5, height = 5)
