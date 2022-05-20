function Kfold_SVM(PH,L,cij,model,K,run)


DT=.1;    % Start-of-event threshold
RT=.8;    % End-of-event threshold

filename=sprintf('%s_%d_fold_%s.mat',model,K,run);


%% Train and validate


tvMCCResults=zeros(7,length(PH)*length(L));
PartResults=zeros(K+1,length(L));
PartResults(1,:)=L;

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
                
                tempMCC=zeros(size(cij));
                for kk=1:length(cij)
                    
               
                    %mdl = fitcsvm(tData,tLabels,'Cost',[0 cij(kk);1 0],'KernelFunction','rbf');
                    mdl = fitcsvm(tData,tLabels,'Cost',[0 cij(kk);1 0]);
                    predictions=predict(mdl,tData);
                    [TP,FP,TN,FN]=getConfMatrix(predictions,tLabels);
                    tempMCC(kk) = (TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
                end
                [~,cidx]=max(tempMCC);
                %mdl = fitcsvm(tData,tLabels,'Cost',[0 cij(cidx);1 0],'KernelFunction','rbf');
                mdl = fitcsvm(tData,tLabels,'Cost',[0 cij(cidx);1 0]);
                predictions=predict(mdl,vData);
                    [TP,FP,TN,FN]=getConfMatrix(predictions,vLabels);
                    PartResults(ll+1,jj)=(TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
                    MCC = MCC+PartResults(ll+1,jj);
                

            end
            MCC=MCC/K;

            tvMCCResults(:,count)=[PH(ii);L(jj);1e-16;K;1e-16;MCC;MCC];
            count=count+1;
            save(filename,'cij','tvMCCResults','PartResults');
        end
        
     
    end
end