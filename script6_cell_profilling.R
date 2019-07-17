suppressPackageStartupMessages({
  suppressWarnings({
    library(tidyverse)
    library(ggpubr)
    library(heatmaply)
  })
})

# set global ggplot theme
theme_set(theme_pubr())

# load data and produce tidy dataframe/tibble

print(str_glue("Reading in and combining features from all cells..."))

imgPaths = paste0(getwd(), "/modelImages/",
                  list.files("modelImages/", pattern = "\\.tif$"))
filesdf = data.frame(img = paste("Img", 1:length(imgPaths)), 
                     paths = imgPaths, stringsAsFactors = F)

filelist = apply(filesdf, 1, function(x) {
  list("img" = x[1], "path" = x[2])
})

excluded = NULL

df = map(filelist, function(img) {
  if (!(img$img %in% excluded)) {
    path = str_glue("modelImages/preprocessed/{img$img}/features_labeled.tsv")
    df = read.table(path, dec = ",")
    df$Image = factor(img$img)
    if (sum(is.na(df$Type)) > 0) {
      print(paste0("Got NAs with:", img$img))
    }
    df = df %>% select(Type, Image, 2:(ncol(df) - 1))
    return(df)
  }
}) %>% plyr::ldply(rbind) %>%
  mutate(Cell.id = paste0("cell.", 1:nrow(.))) %>%
  gather("feature", "value", -c("Type", "Image", 'Cell.id')) %>%
  filter(!(Type %in% c("artifact", "unidentified"))) %>%
  mutate(feature.type = ifelse(startsWith(feature,fixed("x")), 
                               str_extract(feature, 'm|s|b'),
                               str_extract(feature, 'h'))) %>%
  mutate(feature.type = c("m" = "Moment",
                          "s" = "Shape",
                          "b" = "Basic",
                          "h" = "Haralick")[feature.type],
         cell.type = factor(as.character(Type), 
                            levels = c('regular', 'spiky', 'bald', 'ghost'))) %>%
  select(-Type)


my.comparisons = list(
  c("regular", "spiky"),
  c("spiky", "bald"),
  c("bald", "ghost")
)

my_theme = list(
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(face = 'bold'),
        strip.text = element_text(face = "bold", size = 12),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(face = "bold", size= 12))
)
# shape features
shape.features = c("x.0.s.area" = "Area", 
                   "x.0.s.radius.mean" = "Mean Radius")
df %>%
  filter(feature %in% names(shape.features)) %>%
  mutate(feature = factor(shape.features[as.character(feature)])) %>%
  group_by(feature) %>%
  mutate(value = 100*value/max(value)) %>%
  ungroup %>%
  ggplot(aes(cell.type, value, color = cell.type)) +
  facet_wrap(~feature) +
  geom_jitter(alpha=.5) +
  geom_boxplot(color = 'gray20', alpha = 0) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE,
                     comparisons = my.comparisons,tip.length = .01,
                     position = position_nudge(y=-.5),label = "p.signif") +
  stat_compare_means(position = position_nudge(y = .75)) +
  scale_y_log10(breaks = c(5, 10, 25, 50, 100)) +
  coord_cartesian(ylim = c(2.5,1000)) +
  labs(y = "% of maximum value", x = NULL, color = NULL) +
  my_theme


# haralick features
halarick.features = c("h.sen" = "Sen", 
                      "h.ent" = "Ent",
                      "h.den" = "Den")
my.comparisons2 = my.comparisons
my.comparisons2[[3]] = c("regular", "ghost")
df %>%
  filter(feature.type=="Haralick") %>%
  mutate(feature.prefix = str_extract(feature, "h\\.\\w+")) %>%
  group_by(feature.prefix, Cell.id) %>%
  mutate(value = mean(value)) %>%
  ungroup %>%
  filter(feature.prefix %in% names(halarick.features)) %>%
  mutate(feature = factor(halarick.features[as.character(feature.prefix)])) %>%
  group_by(feature) %>%
  mutate(value = 100*value/max(value)) %>%
  ungroup %>%
  ggplot(aes(cell.type, value, color = cell.type)) +
  facet_wrap(~feature) +
  geom_jitter(alpha=.5) +
  geom_boxplot(color='gray20', alpha = 0) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE,
                     comparisons = my.comparisons2,
                     tip.length = .01, label = "p.signif",
                     position = position_nudge(y=-.5),
                     show.legend = FALSE) +
  stat_compare_means(position = position_nudge(y = 25),
                     show.legend = FALSE) +
  # scale_y_log10(breaks = c(50, 60, 70, 80, 90, 100)) +
  scale_y_continuous(breaks = c(50, 60, 70, 80, 90, 100)) +
  coord_cartesian(ylim = c(50,130)) +
  labs(y = "% of maximum value", x = NULL, color = NULL) +
  guides(color = guide_legend(override.aes = list(alpha = 1,
                                                  size = 2))) +
  my_theme

# heatmap with z-scores from wuantile regression or something 
