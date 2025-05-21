feats = readmatrix('feats.csv');
feats = feats(2:11,2);

modelKmeans = load("KNNKmeans.mat").KNNKmeans;
modelSingle = load("KNNSingle.mat").KNNSingle;
modelAverage = load("KNNAverage.mat").KNNAverage;
modelComplete = load("KNNComplete.mat").KNNComplete;
modelWard = load("KNNWard.mat").KNNWard;
modelGMM = load("KNNGMM.mat").KNNGMM;
modelSpect = load("KNNSpect.mat").KNNSpect;
modelHDBSCAN = load("KNNHDBSCAN.mat").KNNHDBSCAN;

IndexPredictor(feats,modelKmeans) %Input features, and model for selected clustering algorithm
