function [TP,FP,TN,FN]=TestThresh(t,d,WT,PN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Input
% t : Time in kyr
% d : Dipole
% WT: Warning threshold (as % of mean dipole magnitude)
%
%      Output
% Count of total number of true positives (TP), false positives (FP), true
% negatives (TN) and false negatives (FN).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Prediction=Forecast(t,d,WT,PN);

TP=0;
FP=0;
TN=0;
FN=0;

for ii=1:length(t)
    if PN(ii)~=2
        if PN(ii)==1 && Prediction(ii)==1
            TP=TP+1;
        elseif PN(ii)==1 && Prediction(ii)==0
            FN=FN+1;
        elseif PN(ii)==0 && Prediction(ii)==1
            FP=FP+1;
        elseif PN(ii)==0 && Prediction(ii)==0
            TN=TN+1;
        end
    end
end
            
        