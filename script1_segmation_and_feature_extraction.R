suppressPackageStartupMessages({
  suppressWarnings({
    library(caTools)
    library(doParallel)
    library(EBImage)
    library(tidyverse)
  })
})

imgPaths = paste0(getwd(), "/modelImages/",
                  list.files("modelImages/", pattern = "\\.tif$"))
filesdf = data.frame(img = paste("Img", 1:length(imgPaths)), 
                     paths = imgPaths, stringsAsFactors = F)

filelist = apply(filesdf, 1, function(x) {
  list("img" = x[1], "path" = x[2])
})

saveRDS(filelist, "filelist.RDS")
cores = 3
cl = makeCluster(cores)
registerDoParallel(cl)
print(str_glue("Image precessing, cell segmentation, and feature extraction for {length(imgPaths)} input files."))

big.list = foreach(i = filelist, 
                   .packages = loadedNamespaces()) %dopar% {
                     img_path = i$path; img_numb = i$img
                     print(str_glue("Processing {img_numb} at \n {img_path}"))
                     img <- readImage(img_path, type = "tiff")
                     img = img[,0:724,]
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
                     feats = computeFeatures(cells[,,1], img[,,1], 
                                             methods.noref = c("computeFeatures.moment", 
                                                               "computeFeatures.shape"),
                                             methods.ref = c("computeFeatures.basic", 
                                                             "computeFeatures.moment")) %>%
                       as.data.frame()
                     haralick = computeFeatures.haralick(cells[,,1], img[,,1])
                     haralick = haralick[rowSums(haralick) > 0,]
                     feats = cbind(feats, haralick)
                     rownames(feats) = paste0("C", 1:nrow(feats))
                     feats$Type <- NA
                     feats = feats[,c(ncol(feats), 1:(ncol(feats) - 1))]
                     # Create results dir
                     if (!dir.exists("modelImages/preprocessed")) {
                       dir.create("modelImages/preprocessed")
                     }
                     processed_files_path = str_glue("modelImages/preprocessed/{img_numb}")
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

saveRDS(big.list, "big.list.RDS")

print(str_glue("Done!"))

