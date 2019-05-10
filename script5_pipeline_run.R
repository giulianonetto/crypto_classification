#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  suppressWarnings({
    library(caret)
    library(tidyverse)
    library(ggpubr)
    library(reshape2)
    library(EBImage)
    library(argparse)
    library(doParallel)
  })
})

vr_help_msg = "Vertical range of pixels to cut image into, once you may have annotations on the bottom of the image. Defaults to 1:724"

# define Argument's parser

parser <- ArgumentParser()

parser$add_argument("-d", "--images_directory", type='character',
                    help="Path to directory where all input .tif images are stored")
parser$add_argument("-c", "--classifier", type='character',
                    help="Path to classifier as a .RDS file.")
parser$add_argument("-o", "--output_directory", type='character',
                    default="output_files",
                    help="Path to directory where all outputs will be stored")
parser$add_argument("-vr", "--vertical_range", type="character",
                    default='1:724',
                    help=vr_help_msg)
parser$add_argument("-t", "--threads", type='integer', default=3,
                    help="Number of threads to use. Defaults to 3.")

args <- parser$parse_args()

images_directory <- args$images_directory
output_directory <- args$output_directory
classifier <- args$classifier
threads <- args$threads
vertical_range <- args$vertical_range


vertical_range = as.integer(str_split(vertical_range, ":", simplify = TRUE)[1,])
vertical_range = seq(vertical_range[1], vertical_range[2])

if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

imgPaths = list.files(images_directory, pattern = "\\.tif$", full.names = T)
filesdf = data.frame(img = paste("Img", 1:length(imgPaths)), 
                     paths = imgPaths, stringsAsFactors = F)

filelist = split.data.frame(filesdf, f = filesdf$img)
filelist = map(filelist, function(x) {
  x[["out"]] = output_directory
  return(x)
})

# set global ggplot theme
theme_set(theme_pubclean(base_size = 16))

# load data
cl = makeCluster(threads)
registerDoParallel(cl)
print(str_glue("Image precessing, cell segmentation, and feature extraction for {length(imgPaths)} input files."))

big.list = foreach(i = filelist, 
                   .packages = loadedNamespaces()) %dopar% {
                     img_path = i$path; img_numb = i$img
                     output_directory = i$out
                     print(str_glue("Processing {img_numb} at \n {img_path}"))
                     img <- readImage(img_path, type = "tiff")
                     if (numberOfFrames(img) == 1) {
                       img = img[,vertical_range]
                     } else {
                       img = img[,vertical_range, 1]
                     }
                     
                     colorMode(img) <- "Grayscale"
                     
                     # Process image for segmentation
                     pix = 651
                     disc = makeBrush(pix, shape = "box")
                     disc = disc / sum(disc)
                     localBackground = filter2(img, disc)
                     # offset = 0.0005
                     offset = 0.01
                     img.binary = (img - localBackground > offset)
                     nucSeed = bwlabel(img.binary)
                     foregroudMask = img - filter2(img, disc) > 0.0003
                     foregroudMask = fillHull(foregroudMask)
                     cells = propagate(img, nucSeed, 
                                       mask = foregroudMask)
                     # filter by object SIZE!
                     min_px_by_obj = 1000
                     tabs = table(cells)
                     toremove = names(tabs[tabs < min_px_by_obj])
                     cells[cells %in% toremove] = 0
                     # compute Features!!
                     feats = computeFeatures(cells, img, 
                                             methods.noref = c("computeFeatures.moment", 
                                                               "computeFeatures.shape"),
                                             methods.ref = c("computeFeatures.basic", 
                                                             "computeFeatures.moment")) %>%
                       as.data.frame()
                     haralick = computeFeatures.haralick(cells, img)
                     haralick = haralick[rowSums(haralick) > 0,]
                     feats = cbind(feats, haralick)
                     rownames(feats) = paste0("C", 1:nrow(feats))
                     # feats$Type <- NA
                     feats = feats[,c(ncol(feats), 1:(ncol(feats) - 1))]
                     # Create results dir
                     results_path = str_glue("{output_directory}/preprocessed")
                     if (!dir.exists(results_path)) {
                       dir.create(results_path)
                     }
                     processed_files_path = str_glue("{results_path}/{img_numb}")
                     dir.create(processed_files_path)
                     # Generate labelled image
                     {
                       png(str_glue("{processed_files_path}/labels.png"),
                           res = 300, width = 1024, height = 724)
                       EBImage::display(colorLabels(cells), method = "raster")
                       text(x = feats$x.0.m.cx, y = feats$x.0.m.cy, 
                            labels = rownames(feats), col = "white", cex = 0.7)
                       dev.off()    
                     }
                     
                     {
                       png(str_glue("{processed_files_path}/labels_painted_objects.png"),
                           res = 300, width = 1024, height = 724)
                       display(paintObjects(cells,
                                            EBImage::toRGB(img), thick = T,
                                            col = "green4"),
                               method = "raster",
                               margin = c(50,50))
                       segments(x0 = feats$x.0.m.cx, x1 = feats$x.0.m.cx,
                                y0 = feats$x.0.m.cy, y1 = feats$x.0.m.cy - 29, 
                                lwd = 1.5,
                                col = "gray15")
                       text(x = feats$x.0.m.cx, y = feats$x.0.m.cy - 40,
                            labels = rownames(feats),
                            col = "gray80",
                            cex = 0.42, font = 2)
                       dev.off()
                     }
                     
                     # Save original cells object and features table
                     saveRDS(cells, str_glue("{processed_files_path}/cells.RDS"))
                     write.table(feats,
                                 str_glue("{processed_files_path}/features.tsv"),
                                 sep = "\t", quote = F, row.names = T,
                                 col.names = T)
                     # Copy original image
                     file.copy(img_path, processed_files_path)
                     if (is.null(cells)) {
                       msg = str_glue("{img_num} returned NULL for cells object.")
                       stop(msg)
                     } else if (is.null(feats)) {
                       msg = str_glue("{img_num} returned NULL for feats object.")
                       stop(msg)
                     } else if (is.null(img)) {
                       msg = str_glue("{img_num} returned NULL for img object.")
                       stop(msg)
                     }
                     
                     # return objects
                     return(list("cells" = cells, "img" = img, "feats" = feats))
                   }

# read in features table
feats = read.table(str_glue("{output_directory}/preprocessed/Img 1/features.tsv"),
                   sep = "\t", header = TRUE)
positions = feats %>% 
  select(x.0.m.cx, x.0.m.cy) %>%
  dplyr::rename("x" = x.0.m.cx, "y" = x.0.m.cy)

feats_to_keep = computeFeatures(big.list[[1]]$cells,
                                big.list[[1]]$img, 
                                properties = TRUE) %>%
  filter(translation.invariant, rotation.invariant) %>%
  select(name) %>% unlist() %>% unname()
df = feats[, colnames(feats) %in% feats_to_keep]

# make predictions
fit = readRDS(classifier)
preds = predict(fit, feats)


# Draw predictions on original images (testset-derived cells only)

print(str_glue("Drawing predictions on original images (testset-derived cells only)..."))


locations = data.frame(x = positions$x, y = positions$y,
                       Cell = rownames(positions),
                       Image = "Img 1",
                       predictions = preds,
                       stringsAsFactors = F)


for (i in filelist) {
  img_id = i$img
  img_path = i$path
  img_numb = as.numeric(str_extract(img_id, "\\d+"))
  subset_loc = locations %>% filter(Image == img_id)
  if (nrow(subset_loc) == 0) next
  files_path = str_glue("{output_directory}/predictions")
  if (!dir.exists(files_path)) {
    dir.create(files_path)
  }
  
  # Generate labelled image
  text_labels = as.character(subset_loc$predictions)
  {
    png(str_glue("{files_path}/Img_{img_numb}_predictions.png"),
        res = 300, width = 1024, height = 724)
    display(paintObjects(big.list[[img_numb]]$cells,
                         EBImage::toRGB(big.list[[img_numb]]$img), 
                         thick = T,
                         col = "green4"),
            method = "raster",
            margin = c(50,50))
    segments(x0 = subset_loc$x, x1 = subset_loc$x,
             y0 = subset_loc$y - 5, y1 = subset_loc$y - 29, 
             lwd = 1.5,
             col = "gray80")
    text(x = subset_loc$x, y = subset_loc$y - 40,
         labels = text_labels,
         col = "gray80",
         cex = 0.42, font = 2)
    dev.off()
  }
}

print(str_glue("Done!"))
