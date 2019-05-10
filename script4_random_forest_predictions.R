suppressPackageStartupMessages({
  suppressWarnings({
    library(caret)
    library(tidyverse)
    library(ggpubr)
    library(reshape2)
    library(EBImage)
  })
})
# define functions
plot_performance <- function(confmat, 
                             classSummary,
                             h=1, v=1,
                             to_include=NULL) {
  if (is.null(to_include)) {
    to_include = c("Sensitivity", "Specificity", "Balanced Accuracy")
    names(to_include) = to_include  
  } else {
    names(to_include) = to_include
  }
  accuracy = round(confmat$overall[["Accuracy"]] * 100)
  mean_params = map(to_include, function(param) {
    x = str_replace(param, " ", "_")
    y = str_glue("Mean_{x}")
    y = round(classSummary[[y]] * 100, 1)
    lab = str_glue("{param}\n(mean = {y} %)")
    return(lab)
  })
  mylables = list("variable" = unlist(mean_params))
  
  p = confmat$byClass %>% as.data.frame() %>%
    dplyr::select(to_include) %>%
    mutate(Type = str_extract(rownames(.), "[a-z]+$")) %>%
    melt() %>% mutate(value = round(value * 100)) %>%
    ggdotchart("Type", "value", fill = "Type",
               color = "Type", add = "segments",
               facet.by = "variable", 
               label = "value",
               panel.labs = mylables,
               font.label = list(color = "white", size = 10,
                                 face = "bold",
                                 vjust = 0.5),
               dot.size = 10, 
               add.params = list("size" = 2.5, 
                                 "color" = "gray40")) + 
    scale_y_continuous(breaks = c(0, 25, 50, 70, 80, 90, 100), 
                       limits = c(0, 105)) +
    labs(x = NULL, y = "Test Performance (%)") +
    ggtitle(str_glue("Overall Accuracy: {accuracy} %")) +
    theme_pubclean() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          strip.text = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = h,
                                     vjust = v, face = "bold",
                                     size = 11)) +
    guides(label = guide_legend(override.aes=list(label=NA)))
  return(p)  
}

# set global ggplot theme
theme_set(theme_pubclean(base_size = 16))

# load data
testset = read.csv("modelImages/preprocessed/testing_set.tsv", 
                   sep = "\t", header = TRUE)
trainset = read.csv("modelImages/preprocessed/training_set.tsv", 
                    sep = "\t", header = TRUE)
big.list = readRDS("big.list.RDS")
fit = readRDS("tuning_fit_random_forest.RDS") # it has all caret's objects

print(str_glue("Making class predictions..."))

# make predictions
preds = predict(fit, testset[, colnames(testset) != "Type"])
conf_mat_rf = confusionMatrix(preds, testset$Type)

# Plot final metrics

print(str_glue("Building performance plot..."))

csum = data.frame(pred = preds,
                  obs = testset$Type) %>% 
  multiClassSummary(lev = levels(testset$Type)) %>% 
  as.list()

conf_plot = plot_performance(confmat = conf_mat_rf, 
                             classSummary = csum) 
plot_path = "modelImages/plots/fig3_predictive_performance.png"
ggsave(plot_path, conf_plot, width = 12, height = 5)


# Confusion Matrix heatmap

print(str_glue("Building Confusion Matrix heatmaps..."))


ordered_levels = c("unidentified", "spiky",
                   "regular", "bald", "artifact")
confusion_dataframe = as.data.frame(conf_mat_rf$table) %>%
  mutate(Reference = factor(Reference, 
                            levels = ordered_levels),
         Prediction = factor(Prediction, 
                             levels = rev(ordered_levels))) %>%
  group_by(Reference) %>%
  mutate(Proportions = round(Freq / sum(Freq),2)*100)

cmat_heatmap_absolute = confusion_dataframe %>%
  ggplot(aes(Prediction, Reference)) + 
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = sprintf("%1.0f", Freq)),
            color = "white", vjust = 0.5) +
  scale_fill_gradient(low = "#00275D", high = "#D70000",
                      limits = c(0, 65),
                      breaks = seq(0, 60, 20)) +
  labs(fill = "Numbers\nof cells") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title = element_text(size = 18, face = "bold", vjust = 0.5),
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.key.width = unit(1.2, "cm"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom") +
  scale_x_discrete(position = "top") +
  ggtitle("Absolute\n")

cmat_heatmap_proportions = confusion_dataframe %>%
  ggplot(aes(Prediction, Reference)) + 
  geom_tile(aes(fill = Proportions)) + 
  geom_text(aes(label = Proportions),
            color = "white", vjust = 0.5) +
  scale_fill_gradient(low = "#00275D", high = "#D70000",
                      breaks = c(0, 25, 50, 75, 100), 
                      limits = c(0, 105)) +
  labs(fill = "Proportions\n(row-wise, %)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title = element_text(size = 18, face = "bold", vjust = 0.5),
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.key.width = unit(1.2, "cm"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom") +
  scale_x_discrete(position = "top") +
  ggtitle("Proportional\n")

cmats_arranged = ggarrange(cmat_heatmap_absolute, 
                           cmat_heatmap_proportions, 
                           ncol = 2)
plots_dir = "modelImages/plots"
ggsave(str_glue("{plots_dir}/fig4_confusion_matrices.png"),
       cmats_arranged, width = 14, height = 6)
ggsave(str_glue("{plots_dir}/fig4a_confusion_matrix_absolute.png"),
       cmat_heatmap_absolute, width = 7.5, height = 6)
ggsave(str_glue("{plots_dir}/fig4b_confusion_matrix_proportional.png"),
       cmat_heatmap_proportions, width = 7.5, height = 6)


# Draw predictions on original images (testset-derived cells only)

print(str_glue("Drawing predictions on original images (testset-derived cells only)..."))

filelist = readRDS("filelist.RDS")
df_bkp = readRDS("df_bkp.RDS")
test_cell_ids = readRDS("cell_ids.RDS")[["test_ids"]]
hits = preds == testset$Type

locations = data.frame(x = df_bkp$x.0.m.cx, y = df_bkp$x.0.m.cy, 
                       Cell = df_bkp$Cell, 
                       Type = factor(df_bkp$Type, 
                                     levels = levels(testset$Type)),
                       Image = df_bkp$Image,
                       stringsAsFactors = F) %>%
  filter(Cell %in% test_cell_ids) %>%
  mutate(predictions = preds, hits = hits)

# example
hit = "gray80"
error = "#FF0A0A"

for (i in filelist) {
  img_id = i$img
  img_path = i$path
  img_numb = as.numeric(str_extract(img_id, "\\d+"))
  subset_loc = locations %>% filter(Image == img_id)
  if (nrow(subset_loc) == 0) next
  files_path = "modelImages/predictions"
  if (!dir.exists(files_path)) {
    dir.create(files_path)
  }
  
  # Generate labelled image
  text_labels = ifelse(subset_loc$hits,
                       as.character(subset_loc$predictions),
                       str_glue("{subset_loc$predictions} ({subset_loc$Type})"))
  {
    png(str_glue("{files_path}/Img_{img_numb}_predictions.png"),
        res = 300, width = 1024, height = 724)
    display(paintObjects(big.list[[img_numb]]$cells[,,1],
                         EBImage::toRGB(big.list[[img_numb]]$img[,,1]), thick = T,
                         col = "green4"),
            method = "raster",
            margin = c(50,50))
    segments(x0 = subset_loc$x, x1 = subset_loc$x,
             y0 = subset_loc$y - 5, y1 = subset_loc$y - 29, lwd = 1.5,
             col = ifelse(subset_loc$hits, hit, error))
    text(x = subset_loc$x, y = subset_loc$y - 40,
         labels = text_labels,
         col = ifelse(subset_loc$hits, hit, error),
         cex = 0.42, font = 2)
    dev.off()
  }
}

print(str_glue("Done!"))
