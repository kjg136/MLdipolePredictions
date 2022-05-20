function [Data,Labels]=GetLabeledData(t,D,DT,RT,PH,WindowSize)

[start,stop]=FindEvents(t,D,DT,RT);  
PN=PosNeg(t,start,stop,PH);

NumberOfWindows=length(PN)-WindowSize+1;
Data=zeros(WindowSize,NumberOfWindows);
Labels=zeros(1,NumberOfWindows);

for ii=1:NumberOfWindows
    Data(:,ii)=D(1+ii-1:WindowSize+ii-1);
    
    if ismember(2,PN(1+ii-1:WindowSize+ii-1))
        Labels(ii)=2;
    elseif ismember(1,PN(1+ii-1:WindowSize+ii-1))
        Labels(ii)=1;
    end
end

Data=transpose(Data);
Labels=transpose(Labels);