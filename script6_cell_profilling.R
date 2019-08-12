suppressPackageStartupMessages({
  suppressWarnings({
    library(tidyverse)
    library(ggpubr)
    library(heatmaply)
    library(pheatmap)
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

# shape features
.features = c(
  "x.0.s.radius.mean" = "Radius Mean", # shape
  "h.con.s1" = "Contrast", "h.con.s2" = "Contrast", # texture
  "x.0.m.majoraxis" = "Major Axis",  # moment
  "x.a.m.majoraxis" = "Major Axis",  # moment
  "x.Ba.m.majoraxis" = "Major Axis", # moment
  "x.a.b.mean" = "Mean Pixel Intensity" # basic
)

dat <- df %>%
  filter(feature %in% names(.features)) %>%
  group_by(Cell.id, feature.type) %>%
  mutate(.value = mean(value),
         .feature = .features[feature]) %>%
  ungroup %>% 
  select(Cell.id, cell.type, .value, .feature) %>%
  unique()

make_plot <- function(.dat, .pars) {
  .dat %>%
    ggplot(aes(cell.type, .value)) +
    geom_jitter(aes(color = cell.type), alpha = .75,
                height = 0, width = .35) +
    stat_summary(geom = 'line', fun.y = mean,
                 color = 'red', size = 1.25,
                 aes(group = .feature)) +
    .pars
}

.features = unique(.features)
names(.features) = .features
annot.positions = c(
  "Major Axis" = 250,
  "Contrast" = 80,
  "Mean Pixel Intensity" = 0.8,
  "Radius Mean" = 100
)
plt.list <- map(.features, function(.feat) {
  df = dat %>% filter(.feature == .feat)
  .ylims = c(min(df$.value)*.5, max(df$.value)*1.5)
  .pars = list(
    stat_compare_means(label.y = annot.positions[.feat]),
    stat_compare_means(comparisons = my.comparisons,
                       tip.length = .02,
                       label = 'p.signif'),
    ylim(.ylims),
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = 'bold'),
          plot.title = element_text(face = "bold", size = 14),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(face = "bold", size= 12)),
    labs(x=NULL, y=NULL, color=NULL),
    ggtitle(.feat)
  )
  p <- make_plot(.dat = df, .pars = .pars)
  return(p)
})

.feat_plot <- ggarrange(plotlist = plt.list, ncol = 2, nrow = 2,
                        legend = 'right', common.legend = TRUE)
ggsave('modelImages/features_plot.png',.feat_plot,
       width = 10, height = 6)
# heatmap with z-scores from wuantile regression or something 

all.features = df$feature %>% unique()
# remove translation and rotation variant features
ix = !str_detect(all.features,"c[xy]|theta") 
all.features = all.features[ix]
names(all.features) = all.features

model.data <- map(all.features, function(.feat) {
  .df = df %>%
    filter(feature == .feat) %>%
    select(Cell.id, feature, value, cell.type) %>%
    spread(feature, value)
  fit = MASS::rlm(rlang::eval_tidy(as.name(.feat)) ~ cell.type,
                  data = .df)
  .coefs = coef(fit)
  pvals = sapply(names(.coefs), function(i) {
    sfsmisc::f.robftest(fit, var = i)$p.value
                        })
  estimates = .coefs
  estimates[2:4] = sapply(estimates[2:4], 
                          function(i) estimates[1] + i )
  
  coefs_to_cells = c("regular", "spiky", "bald", "phantom")
  names(coefs_to_cells) = names(.coefs)
  
  out.df = data.frame(
    "estimates" = estimates,
    ".coefs" = .coefs,
    "p.value" = pvals,
    "cell.type" = coefs_to_cells[names(.coefs)]
  )
  return(out.df)
}) 

model.data <- model.data %>% plyr::ldply(rbind) %>%
  mutate(p.value = p.adjust(p.value, method = "BH")) %>%
  mutate(corrected_estimates = {
    ifelse(p.value <= 0.05, estimates, estimates - .coefs)
    })
  
hm.data = model.data %>%
  filter(.id!="x.Ba.b.q05") %>% # same for all cell types, breaks scaling
  select(.id, cell.type, corrected_estimates) %>%
  spread(cell.type, corrected_estimates)


rownames(hm.data) = hm.data$.id
df = df %>% select(feature, feature.type) %>% unique()
rownames(df) = df$feature
.feat_types = df[hm.data$.id, "feature.type"] %>%
  data.frame(row.names = rownames(hm.data))
hm.data$.id = NULL
h.plot = pheatmap(hm.data, annotation_row = .feat_types, scale = "row")

save_pheatmap_png <- function(x, filename, width=2800, height=2500, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(h.plot, "modelImages/heatmap.png")
