clear;
addpath('polyatreetest');
addpath('incremental-ks-master');

%% parameters
rng(1);     % For reproducibility
windowsize = 100;
SVMkernel = 'gaussian';

testmaxnum = -1;

%fid = fopen('incremental-ks-master/datasets/Bike_num.data');
%fid = fopen('incremental-ks-master/datasets/Arabic_Paper_Fixed.data');
%fid = fopen('incremental-ks-master/datasets/Insects_num.data');
fid = fopen('incremental-ks-master/datasets/Posture_num.data');
%fid = fopen('incremental-ks-master/datasets/Keystroke_norm.data');


%% read top windowsize data points to build initial classifier

% read data

data = [];
datanum = 0;

tline = fgetl(fid);

while ischar(tline)
    %disp(tline)
    tmp     = cell2mat(cellfun(@str2num, strsplit(tline,','),'un',0));
    data    = [data ; tmp];
    datanum = datanum + 1;
    
    %pause    
    if datanum == windowsize
        break;
    else
        tline = fgetl(fid);
    end
end

featurenum = size(data, 2) -1;

% train classifier

referenceX = data(:, 1:featurenum);
referenceY = data(:, end);

%%%%%% SVM  'gaussian', 
svmtemp  = templateSVM('KernelFunction',SVMkernel);
SVMModel = fitcecoc(referenceX, referenceY, 'Learners',svmtemp, 'FitPosterior',1);
insampleerror = resubLoss(SVMModel)

%%%%%% kNN   Distance = {minkowski, euclidean, chebychev}
%KNNModel = fitcknn(referenceX, referenceY, 'Distance','minkowski', 'NumNeighbors',5);
%insampleerror = resubLoss(KNNModel)


%% read data point one by one from windowsize+1 

classification_result = 0;
testdatanum = 0;
driftnum    = 0;

currentX = referenceX;
currentY = referenceY;

tline = fgetl(fid);
while ischar(tline)
    
    % obtain next datum
    %disp(tline)
    testdatanum   = testdatanum + 1    
    newdatum      = cell2mat(cellfun(@str2num, strsplit(tline,','),'un',0));  
    
    % predict new coming datum      
    [label,NegLoss,PBScore,Posterior] = predict(SVMModel,newdatum(1:featurenum));
    %[label,score,cost] = predict(KNNModel,newdatum(1:featurenum));
    
    if label == newdatum(end)    
        classification_result   = classification_result + 1;
        %disp('correct');
    else
        %disp('wrong');
    end
           
    %% detect drift
    currentX(1, :) = [];
    currentY(1)    = [];
    currentX       = [currentX; newdatum(1:featurenum)];
    currentY       = [currentY; newdatum(end)];
    
    drift = 0;  
    f     = 1;
    while ~drift && f <= featurenum
        
        % using KS test
        %[drift, p] = kstest2(referenceX(:, f), currentX(:, f), 'alpha', 0.05);
        
        % using Poly-tree test
        %[drift, post, stats] = PTtest(referenceX(:, f), currentX(:, f), 'partition', 'empirical');
        [drift, post, stats] = PTtest(referenceX(:, f), currentX(:, f), 'normalize', true);

        f = f + 1;        
    end
    
    %drift = 1;
    %%
    if drift
        fprintf('%d drift! \n', testdatanum);
        driftnum   = driftnum + 1;
        
        % update reference set and train a new classifier
        referenceX = currentX;
        referenceY = currentY;
        SVMModel   = fitcecoc(referenceX, referenceY, 'Learners',svmtemp, 'FitPosterior',1);
        %KNNModel   = fitcknn(referenceX, referenceY, 'Distance','minkowski', 'NumNeighbors',5);
    end
    
    % read next datum    
    if testmaxnum == testdatanum
        break;
    else
        tline = fgetl(fid);
    end
end 

fclose(fid);

%% save results 

labelcorrect_rate = classification_result/testdatanum
drift_rate        = driftnum/testdatanum

%save topliner.mat labelcorrect_rate drift_rate;
