suppressPackageStartupMessages({
  suppressWarnings({
    library(tidyverse)
    library(ggfortify)
    library(ggpubr)
    library(EBImage)
    library(Rtsne)
    library(vegan)
  })
})


theme_set(theme_pubr())

print(str_glue("Reading in and combining features from all cells..."))

imgPaths = paste0(getwd(), "/modelImages/",
                  list.files("modelImages/", pattern = "\\.tif$"))
filesdf = data.frame(img = paste("Img", 1:length(imgPaths)), 
                     paths = imgPaths, stringsAsFactors = F)

filelist = apply(filesdf, 1, function(x) {
  list("img" = x[1], "path" = x[2])
})

big.list = readRDS("big.list.RDS")
excluded = c("none")
# read in features tables
feats_list = lapply(filelist, function(img) {
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
  
})

saveRDS(feats_list, "feats_list.RDS")

df = plyr::ldply(feats_list, rbind)
df$Cell = paste0("Cell", rownames(df))
df_bkp = df
saveRDS(df_bkp, "df_bkp.RDS")
df = df[, -2]
props =  EBImage::computeFeatures(big.list[[1]]$cells[,,1], 
                                  big.list[[1]]$img[,,1],
                                  properties = TRUE)
to.remove <- props[!(props$translation.invariant & props$rotation.invariant),
                   "name"]
df <- df %>% 
  select(-to.remove) %>%
  filter(Type != "ghost") %>%
  mutate(Type = factor(as.character(Type),
                       levels = c("bald",
                                  "spiky",
                                  "regular",
                                  "artifact",
                                  "unidentified")))
df$Type = droplevels(df$Type, "ghost")


# Store table with features from all labeled cells / images
all_features_path = "features_all_labeled.tsv"
print(str_glue("Saved all features at: \n {all_features_path}"))
write.table(df,
            all_features_path,
            sep = "\t", quote = F, row.names = T,
            col.names = T)

print(str_glue("Plotting Principal Component Analysis..."))


pca = prcomp(df[, -c(1,ncol(df))], scale = T)

p1 = ggplot2::autoplot(pca, data = df, 
                       colour = "Type", size = 3) +
  theme(plot.title = element_text(hjust = 0.475, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 18, face = 'bold'),
        legend.title = element_text(hjust = 0, 
                                    face = "bold")) +
  labs(color = NULL) +
  ggtitle("Front View") +
  scale_color_brewer(type = "qual",palette = 6)

p2 = ggplot2::autoplot(pca,x = 1, y = 3, data = df,
                       colour = "Type", size = 3) +
  theme(plot.title = element_text(hjust = 0.475, face = "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 18, face = 'bold'),
        legend.title = element_text(hjust = 0, face = "bold")) +
  labs(color = NULL) +
  ggtitle("Top View") +
  scale_color_brewer(type = "qual",palette = 6)

pcaplot = ggarrange(p1, p2, ncol = 2, 
                    common.legend = T, legend = "right")
plots_dir = "modelImages/plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

pca_plot_path = str_glue("{plots_dir}/fig1_pca.jpeg")
print(str_glue("Saved PCA plot at: \n {pca_plot_path}"))
ggsave(pca_plot_path, pcaplot, width = 16, height = 5.5)

# Run below for interactive, 3D visualization
#plotly::plot_ly(data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3, color = df$Type)

# PERMANOVA only with labeled cells
cell_types = c("bald", "regular", "spiky")
distance_matrix = dist(df[df$Type %in% cell_types, -c(1,2, ncol(df))])
metadata = data.frame("Type" = df[df$Type %in% cell_types, "Type"])
permanova = vegan::adonis(distance_matrix ~ Type, data = metadata)
saveRDS(permanova, "permanova_test.RDS")

# t-sne

dat = df[, -c(1,ncol(df))]

set.seed(9)  
tsne = Rtsne(as.matrix(dat), check_duplicates = FALSE, 
             max_iter = 1000,
             pca = TRUE, perplexity = 200, theta=.5, 
             dims=3, num_threads = 0)

## getting the two dimension matrix
tsne.axis = as.data.frame(tsne$Y)
tsne.axis$type = df$Type
p1p2=ggplot(tsne.axis, aes(V1, V2, color = type)) +
  geom_point(size=3, alpha=.4)+
  stat_ellipse(aes(group=type),level = .5, size =2)
p1p3=ggplot(tsne.axis, aes(V1, V3, color = type)) +
  geom_point(size=3)+
  stat_ellipse(level = .5, size = 2)
ggarrange(p1p2,p1p3, ncol=2, common.legend = T)

# mds

dist.dat = dist(dat, method = "minkowski")
mds = cmdscale(dist.dat, k = 3)
mds.dat = data.frame(mds, type = df$Type)

p1p2=ggplot(mds.dat, aes(X1, X2, color = type)) +
  geom_point(size=3)
p1p3=ggplot(mds.dat, aes(X1, X3, color = type)) +
  geom_point(size=3)
ggarrange(p1p2,p1p3, ncol=2, common.legend = T)




