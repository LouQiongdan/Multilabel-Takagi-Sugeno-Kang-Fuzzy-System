clear all
clc
tic
clear

load('./MLTSKFS/data/CAL500.mat');

target(find(target==-1))=0;
data = (mapminmax(data',0,1))';
oldOptmParameter=struct('alpha_searchrange',[0.01,0.1,1,10,100],'beta_searchrange',[0.01,0.1,1,10,100],'gamma_searchrange',[0.01,0.1,1,10,100],...
    'maxIter',100,'minimumLossMargin',0.01,'outputtempresult',0,'drawConvergence',0,'bQuiet',0);
TSKoptions=struct('k_searchrange',[2,3],'h_searchrange',[0.1,1,10,100]);

[ BestParameter, BestResult ] = ML_TSKFS_adaptive_validate( data, target, oldOptmParameter,TSKoptions);

toc 



