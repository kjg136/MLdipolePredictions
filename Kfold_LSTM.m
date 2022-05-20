function Kfold_LSTM(PH,L,cij,model,K,run)

hsize=100;
Nepoch=10;
DT=.1;    % Start-of-event threshold
RT=.8;    % End-of-event threshold

filename=sprintf('%s_%d_fold_%s.mat',model,K,run);

options = trainingOptions('adam', ...
    'MaxEpochs',Nepoch, ...
    ...
    'InitialLearnRate', 0.01, ...
    'GradientThreshold', 1, ...
    'ExecutionEnvironment',"auto",...
    ...'plots','training-progress', ...
    'Verbose',false);
%'MiniBatchSize', 150, ...
%'SequenceLength', 1000, ...
% 'plots','training-progress', ...

catnames=["N","Y"];


%% Train and validate


tvMCCResults=zeros(7,length(PH)*length(L));
PartResults=zeros(K+1,length(L));
PartResults(1,:)=L;
tMCC=zeros(length(cij)+4,length(PH)*length(L)*K);

count=1;
count2=1;

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
                
                tData=num2cell(tData,2);
                    tLabels=categorical(tLabels,[0 1],catnames);
                
                vP=PData(Pidx(ll,1):Pidx(ll,2),:);
                vN=NData(Nidx(ll,1):Nidx(ll,2),:);
                
                vData=[vP;vN];
                vLabels=[ones(size(vP,1),1);zeros(size(vN,1),1)];
                
                 vData=num2cell(vData,2);
                 vLabels=categorical(vLabels,[0 1],catnames);
                
                tempMCC=zeros(size(cij));
                for kk=1:length(cij)
                    
                     layers = [ ...
                    sequenceInputLayer(1)
                    bilstmLayer(hsize,'OutputMode','last')
                    fullyConnectedLayer(2)
                    softmaxLayer
                    classificationLayer('Classes',categorical(catnames),'ClassWeights',[cij(kk) 1])
                ];
                
                
                
                net = trainNetwork(tData,tLabels,layers,options);
                
                predictions = classify(net,vData);
                    [TP,FP,TN,FN]=getConfMatrixLSTM(predictions,vLabels);
                    tempMCC(kk) = (TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
                end
                tMCC(:,count2)=[PH(ii);L(jj);1e-16;K;tempMCC.'];
                count2=count2+1;
                [~,cidx]=max(tempMCC);
                
                layers = [ ...
                    sequenceInputLayer(1)
                    bilstmLayer(hsize,'OutputMode','last')
                    fullyConnectedLayer(2)
                    softmaxLayer
                    classificationLayer('Classes',categorical(catnames),'ClassWeights',[cij(cidx) 1])
                ];
                
            net = trainNetwork(tData,tLabels,layers,options);
                
               predictions = classify(net,vData);
                    [TP,FP,TN,FN]=getConfMatrixLSTM(predictions,vLabels);
                    PartResults(ll+1,jj)=(TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
                    MCC = MCC+PartResults(ll+1,jj);
                

            end
            MCC=MCC/K;

            tvMCCResults(:,count)=[PH(ii);L(jj);1e-16;K;1e-16;MCC;MCC];
            count=count+1;
            save(filename,'cij','tvMCCResults','PartResults','tMCC');
        end
        
     
    end
end