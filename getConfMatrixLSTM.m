function [TP,FP,TN,FN]=getConfMatrix(predictions,Labels)

TP=sum(predictions=="Y" & Labels=="Y");
FP=sum(predictions=="Y" & Labels=="N");
TN=sum(predictions=="N" & Labels=="N");
FN=sum(predictions=="N" & Labels=="Y");
