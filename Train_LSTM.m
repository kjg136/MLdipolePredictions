function Train_LSTM(PH,L,cij,RusRatio,model,TrainingEvents,run)

rng(1);
DT=.1;    % Start-of-event threshold
RT=.8;    % End-of-event threshold

filename=sprintf('%s_%d_%s_LSTM.mat',model,TrainingEvents,run);

if TrainingEvents~=0
    if strcmp(model,'ThreeD') || strcmp(model,'ThreeDCoarse')
        load(model)
    else
        load(strcat(model,'b'));
    end
    [start,~]=FindEvents(t,d,DT,RT);
    idx=find(t==start(TrainingEvents,1));
    tt=t(1:idx);
    tD=d(1:idx);
else
    load(model)
    tt=t;
    tD=d;
end

dt=tt(2)-tt(1);

if TrainingEvents~=0 && (strcmp(model,'ThreeD') || strcmp(model,'ThreeDCoarse'))
    vD=d(idx+1:end);
    vt=t(idx+1:end)-t(idx+1);
else
    load(model)
    vD=d;
    vt=t;
end

%% Train and validate

MCCResults=zeros(5,length(PH));
tvMCCResults=zeros(7,length(PH)*length(L));
FullResults=zeros(4+length(cij),length(PH)*length(L)*length(RusRatio));
count=1;
count2=1;
catnames=["N","Y"];

for ii=1:length(PH)
    for jj=1:length(L)
        jj
        Lscore=0;
        WindowSize=round(L(jj)/dt);
        [Data,Labels]=GetLabeledData(tt,tD,DT,RT,PH(ii),WindowSize);
        Data=abs(Data);
        
        FulltData=num2cell(Data,2);
        %FulltLabels=discretize(Labels,[-.5 .5 1.5],'categorical',catnames);
        FulltLabels=categorical(Labels,[0 1],catnames);
        
        
        
        [vData,vLabels]=GetLabeledData(vt,vD,DT,RT,PH(ii),WindowSize);
        vData=abs(vData);
        vData=num2cell(vData,2);
        
        %vLabels=discretize(vLabels,[-.5 .5 1.5],'categorical',catnames);
        vLabels=categorical(vLabels,[0 1],catnames);
        
        %% Sort and train
        
        NP=sum(Labels==1);
        NN=sum(Labels==0);
        
        for ll=1:length(RusRatio)
            TrainingSizeP=NP;
            TrainingSizeN=min(NN,floor(TrainingSizeP/RusRatio(ll)));
            TrainingSize=TrainingSizeP+TrainingSizeN;
            
            tData=zeros(TrainingSize,WindowSize);
            
            idx=find(Labels==0);
            idx=idx(randperm(length(idx)));
            tData(1:TrainingSizeN,:)=Data(idx(1:TrainingSizeN),:);
            
            idx=find(Labels==1);
            idx=idx(randperm(length(idx)));
            tData(TrainingSizeN+1:end,:)=Data(idx(1:TrainingSizeP),:);
            
            tData=num2cell(tData,2);
            
            tLabels=[zeros(TrainingSizeN,1);ones(TrainingSizeP,1)];
            
            %tLabels=discretize(tLabels,[-.5 .5 1.5],'categorical',catnames);
            tLabels=categorical(tLabels,[0 1],catnames);
            
            FullResults(1:4,count)=[PH(ii);L(jj);RusRatio(ll);TrainingEvents];
            
            for kk=1:length(cij)
                
                layers = [ ...
                    sequenceInputLayer(1)
                    bilstmLayer(100,'OutputMode','last')
                    fullyConnectedLayer(2)
                    softmaxLayer
                    classificationLayer('Classes',categorical(catnames),'ClassWeights',[cij(kk) 1])
                ];
                
                options = trainingOptions('adam', ...
                    'MaxEpochs',10, ...
                    ...
                    'InitialLearnRate', 0.01, ...
                    'GradientThreshold', 1, ...
                    'ExecutionEnvironment',"auto",...
                    ...'plots','training-progress', ...
                    'Verbose',false);
                %'MiniBatchSize', 150, ...
                %'SequenceLength', 1000, ...
                % 'plots','training-progress', ...
                
                net = trainNetwork(tData,tLabels,layers,options);
                
                predictions = classify(net,FulltData);%,'SequenceLength',1000);
                
                [TP,FP,TN,FN]=getConfMatrixLSTM(predictions,FulltLabels);
                MCC = (TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
                
                if MCC>MCCResults(end,ii)
                    MCCResults(:,ii)=[PH(ii);L(jj);RusRatio(ll);cij(kk);MCC];
                end
                
                if MCC>Lscore
                    Lscore=MCC;
                    tvMCCResults(1:end-1,count2)=[PH(ii);L(jj);RusRatio(ll);TrainingEvents;cij(kk);MCC];
                    netOpt=net;
                end
                
                FullResults(4+kk,count)=MCC;
                save(filename,'MCCResults','FullResults','cij','tvMCCResults');
            end
            count=count+1;
        end
        
        if Lscore~=0
            tic
            predictions = classify(netOpt,vData);
            toc
            [TP,FP,TN,FN]=getConfMatrixLSTM(predictions,vLabels);
            MCC = (TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
            tvMCCResults(end,count2)=MCC;
        end
        save(filename,'MCCResults','FullResults','cij','tvMCCResults');
        count2=count2+1;
    end
end

