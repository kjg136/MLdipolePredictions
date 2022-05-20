function MCC=Kfold_Thresh(PH,WarningThresh,model,K,run)

L=1;
DT=.1;    % Start-of-event threshold
RT=.8;    % End-of-event threshold

filename=sprintf('%s_%d_fold_Thresh_%s.mat',model,K,run);


%% Train and validate

MCCResults=zeros(5,length(PH));
tvMCCResults=zeros(7,length(PH)*length(L));
FullResults=zeros(4+length(WarningThresh),length(PH)*length(L));
count=1;

load(model)
dt=t(2)-t(1);

for ii=1:length(PH)
    for jj=1:length(L)
        WindowSize=round(L(jj)/dt);
        [Data,Labels]=GetLabeledData(t,d,DT,RT,PH(ii),WindowSize);
        Data=abs(Data);
        
        %% Sort and train
        
        rng(1);
        NP=sum(Labels==1);
        Pidx=transpose(1:round(NP/K):K*round(NP/K));
        Pidx=[Pidx,[Pidx(2:end)-1;NP]];
        
        NN=sum(Labels==0);
        Nidx=transpose(1:round(NN/K):K*round(NN/K));
        Nidx=[Nidx,[Nidx(2:end)-1;NN]];
        
        idx=find(Labels==0);
        idx=idx(randperm(length(idx)));
        NData=Data(idx,:);
        
        idx=find(Labels==1);
        idx=idx(randperm(length(idx)));
        PData=Data(idx,:);
        
        
        FullResults(1:4,count)=[PH(ii);L(jj);1e-16;K];
        MCC=0;
        for ll=1:K
            tP=PData(setdiff(1:NP,Pidx(ll,1):Pidx(ll,2)),:);
            tN=NData(setdiff(1:NN,Nidx(ll,1):Nidx(ll,2)),:);
            tData=[tP;tN];
            tLabels=[ones(size(tP,1),1);zeros(size(tN,1),1)];
            
            vP=PData(Pidx(ll,1):Pidx(ll,2),:);
            vN=NData(Nidx(ll,1):Nidx(ll,2),:);
            
            vData=[vP;vN];
            vLabels=[ones(size(vP,1),1);zeros(size(vN,1),1)];
            
            tempMCC=zeros(size(WarningThresh));
            for kk=1:length(WarningThresh)
                
                WT=WarningThresh(kk);
                [TP,FP,TN,FN]=TestThresh(tData,tData,WT,tLabels);
                tempMCC(kk)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
            end
            [~,WTidx]=max(tempMCC);
            WT=WarningThresh(WTidx);
            [TP,FP,TN,FN]=TestThresh(vData,vData,WT,vLabels);
            MCC=MCC+(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
            
            %save(filename,'MCCResults','FullResults','WarningThresh','tvMCCResults');
        end
        MCC=MCC/K;
        tvMCCResults(:,count)=[PH(ii);L(jj);1e-16;K;1e-16;MCC;MCC];
        count=count+1;
        %save(filename,'MCCResults','FullResults','WarningThresh','tvMCCResults');
    end
    
    
end
MCC=tvMCCResults(end);

