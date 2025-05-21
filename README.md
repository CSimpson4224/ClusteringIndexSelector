# Clustering Vaidity Index Selector

This repository contains the 8 meta-models for selecting an apropriate clustering validity index based on the selected clustering algorithm, based on the paper "Instance Space of Clustering Validation".
These models predict the performance for the following 9 validity indexes; Wemmert Gancarski, Point-Biserial, AUCC, Silhouette, CDbw, VRC, DB, WB and DBCV.
Models are available for the following clustering algorithms; K-Means, Single-Linkage, Average-Linkage, Complete-Linkage, Ward-Linkage, EM-GMM, Spectral Clustering, and HDBSCAN*. 

The DataFeatures function calculates 10 select meta-features for a given clustering dataset, which are then input into the IndexPredictor function alongside the KNN model file for the selected clustering algorithm.
The output of this is which of the 9 assessed validity indexes are predicted to perform well, along side a selected index based on the most likley correct prediction.

The files MakeFeats and MakePrediction provide an example useage with the provided test dataset.
