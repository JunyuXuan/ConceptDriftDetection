
clear;
addpath('polyatreetest');
addpath('incremental-ks-master');

tic

%% parameters
rng(1);     % For reproducibility
windowsize = 400;
SVMkernel = 'gaussian';

c = 1;
threshold = 0.5;
normalize = true;
estimate_c = false;
partition = [];

testmaxnum = -1;

fid = fopen('incremental-ks-master/datasets/Keystroke.data');%1600

load keystroke_f_m;

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

%% intitial Polya tree

PT = cell(1, featurenum);

for f = 1 : featurenum
    [LOR, PTf] = iniPT(Keystroke_f_m(f, 1), Keystroke_f_m(f, 2), referenceX(:, f), c);
    PT{f}      = PTf; 
end

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
    drift = 0;  
    f     = 1;
    while f <= featurenum
        
        % using KS test
        % [drift, p] = kstest2(referenceX(:, f), currentX(:, f), 'alpha', 0.05);
        
        % using Poly-tree test
        %[drift, post, stats] = PTtest(referenceX(:, f), currentX(:, f), 'partition', 'empirical');
        %[drift, post, stats] = PTtest(referenceX(:, f), currentX(:, f), 'c', c, 'threshold', threshold, 'normalize', normalize, 'estimate_c', estimate_c, 'partition', partition);
        [fdrift, post, stats, PTf] = incPTtest(PT{f}, currentX(1, f), newdatum(f), c, threshold);
        drift = drift + fdrift;
        PT{f} = PTf;
        f     = f + 1;        
    end
    
    currentX(1, :) = [];
    currentY(1)    = [];
    currentX       = [currentX; newdatum(1:featurenum)];
    currentY       = [currentY; newdatum(end)];
            
    %%
    if drift > 0
        fprintf('%d drift! \n', testdatanum);
        driftnum   = driftnum + 1;
        
        % update reference set and train a new classifier
        referenceX = currentX;
        referenceY = currentY;
        SVMModel   = fitcecoc(referenceX, referenceY, 'Learners',svmtemp, 'FitPosterior',1);
        %KNNModel   = fitcknn(referenceX, referenceY, 'Distance','minkowski', 'NumNeighbors',5);
        
        % update PT
        for f = 1 : featurenum            
            [LOR, PTf] = iniPT(Keystroke_f_m(f, 1), Keystroke_f_m(f, 2), referenceX(:, f), c);
            PT{f}      = PTf;
        end
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

toc
%save topliner.mat labelcorrect_rate drift_rate;