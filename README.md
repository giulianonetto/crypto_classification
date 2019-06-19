# Automated machine learning pipeline predicts Cryptococcus gattii morphotypes in Electron Microscopy images

Code and data for paper describing ML-based classification of Cryptococcus gattii cell images.

* __Script 1__ (script1_segmation_and_feature_extraction.R): reads in the image files, processes each file, performs cell segmentation, and extracts many features for each cell (as described [here](https://www.bioconductor.org/packages/devel/bioc/manuals/EBImage/man/EBImage.pdf) for the computeFeatures() function from the EBImage package available from Baicondutor). It takes about one minute using 3 threads (8GB RAM notebook).

* __Script 2__ (script2_combine_features_and_pca.R): takes manually labeled feature tables (e.g. modelImages/preprocessed/Img\ 9/features_labeled.tsv for image 9), combines them together, and plots Principal Componente Analysis (saved in modelImages/plots/fig1_pca.png). 

* __Script 3__ (script3_random_forest_tuning.R): splits data into training and testing sets and performs model tuning (with grid search) for Random Forest classifier building with the [Caret](http://topepo.github.io/caret/index.html) package.

* __Script 4__ (script4_random_forest_predictions.R): using the top-performing model from Script3, it makes class predictions for cells in the test set, builds performance plot (Sensitivity, Specificity, and Balanced Accuracy measurements), plots confusion matrix, and draws the predictions on the original images (testset-derived cells only).
* __Script 5__ (script5_other_images.R): plots supp figs and heatmap of class probabilities.

  The required packages are listed in the __requirements.txt__ file.

