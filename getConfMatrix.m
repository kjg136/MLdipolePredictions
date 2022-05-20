function [TP,FP,TN,FN]=getConfMatrix(predictions,Labels)

TP=sum(predictions==1 & Labels==1);
FP=sum(predictions==1 & Labels==0);
TN=sum(predictions==0 & Labels==0);
FN=sum(predictions==0 & Labels==1);
