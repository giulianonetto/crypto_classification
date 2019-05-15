suppressPackageStartupMessages({
  suppressWarnings({
    library(tidyverse)
    library(ggpubr)
    library(reshape2)
    library(EBImage)
    library(heatmaply)
    library(randomForest)
  })
})

# set global ggplot theme
theme_set(theme_pubr(base_size = 16))

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
# Figure 6
h = heatmaply(df, file = filepath, 
              showticklabels = c(T, F), 
              width = 700, height = 750, fontsize_col = 14)
htmlwidgets::saveWidget(h,"heatmap_probabilities.html")
system("mv heatmap_probabilities.html modelImages/plots")

set.seed(5996)
fit_rf = randomForest(trainset[,-1],
                      trainset$Type,
                      ntree = 2000,
                      nodesize = fit$bestTune$minNode,
                      mtry = fit$bestTune$predFixed,
                      proximity = T,
                      localImp = T,
                      keep.forest = T)

# Fig 7
imp_df = importance(fit_rf, type = 2) %>% 
  data.frame() %>%
  mutate(feature = str_replace(rownames(.), "x.a.", "intensity - ")) %>%
  mutate(feature = str_replace(feature, "x.0.", "binary - ")) %>%
  mutate(feature = str_replace(feature, "x.Ba.", "top hat - "),
         type = str_extract(feature, " \\- b| \\- m|^h| \\- s")) %>%
  mutate(type = factor(type, labels = c("Basic", "Moment", 
                                        "Shape", "Haralick"))) %>%
  arrange(-desc(MeanDecreaseGini)) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

imp_plot = imp_df %>%
  ggplot(aes(feature, MeanDecreaseGini)) +
  geom_col(aes(fill = type), col = "gray20") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(face = "bold")) +
  labs(x = "Features", y = "Mean Decrease Gini",
       fill = NULL) +
  scale_fill_manual(values = c("#00A84C", "#2572A1",
                               "#A12525", "#8C4F02"))
ggsave("modelImages/plots/fig8_variable_importance.jpeg",
       imp_plot, width = 7, height = 9)

# Supp Fig 1 - [proximity plot]
rf_prox = cmdscale(1 - fit_rf$proximity) %>%
  data.frame("type" = trainset$Type)
prox_plot = rf_prox %>%
  ggplot() +
  geom_point(aes(x = X1, y = X2, col = type),
             size = 3, alpha = 0.6) +
  scale_color_brewer(type = "qual", palette = 2) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(col = NULL, x = "Axis1", y = "Axis2")
plotpath = "modelImages/plots/supp_fig1_proximity_plot.jpeg"
ggsave(plotpath, prox_plot, width = 7, height = 5.5)

# Supp Fig2 - top 5 features correlations

top_5_features = rev(levels(imp_df$feature))[1:5] %>%
  str_replace("intensity - ", "x.a.") %>%
  str_replace("binary - ", "x.0.") %>%
  str_replace("top hat - ", "x.Ba.")

top_5_feat_cor_plot = GGally::ggpairs(trainset, 
                                      columns = top_10_features,
                                      mapping = aes(colour = Type,
                                                    alpha = .5)) +
  theme(strip.text = element_text(face = "bold", size = 12))
ggsave("modelImages/plots/supp_fig2_top_5_feats_correlations.jpeg",
       top_5_feat_cor_plot, width = 12, height = 9)
