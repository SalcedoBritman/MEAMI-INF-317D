clear all;
close all;
fclose('all');
clc;

PATHIN = input('Entre com o path do STUDY: ','s');

load([PATHIN filesep 'picoAlfaADDnoAvg.mat']); 
%picoAlfaADD : cell array (16 subjects, 33 epochs x 19 features)
load([PATHIN filesep 'freqAlfaADDnoAvg.mat']); 
%freqAlfaADD : cell array (16 subjects, 33 epochs x 19 features)
load([PATHIN filesep 'picoAlfaNOLDnoAvg.mat']); 
%picoAlfaNOLD : cell array (19 subjects, 33 epochs x 19 features)
load([PATHIN filesep 'freqAlfaNOLDnoAvg.mat']); 
%freqAlfaNOLD : cell array (19 subjects, 33 epochs x 19 features)

%nSubj = length(picoAlfaADD) + length(picoAlfaNOLD) - 3;
nSubj = length(picoAlfaADDnoAvg) + length(picoAlfaNOLDnoAvg) - 3;

allSubjPicos = zeros(size(picoAlfaADDnoAvg{1},2),size(picoAlfaADDnoAvg{1},1),nSubj);
allSubjFreqs = zeros(size(picoAlfaADDnoAvg{1},2),size(picoAlfaADDnoAvg{1},1),nSubj);
allSubjTot = zeros(2*size(picoAlfaADDnoAvg{1},2),size(picoAlfaADDnoAvg{1},1),nSubj);
allSubjLabels = [];
for nf = 1:nSubj
 allSubjLabels = [allSubjLabels blanks(size(picoAlfaADDnoAvg{1},1))']; % epoch labels
end

nSubjVec = 1:nSubj;
group = blanks(nSubj);

for ns = 1:nSubj
    
    if ns < 17 % ADD
      for e = 1:size(picoAlfaADDnoAvg{ns},1) % epochs
        allSubjLabels(e,ns) = 'A';  % all epochs labels
      end
      allSubjPicos(:,:,ns) = picoAlfaADDnoAvg{ns}';
      allSubjFreqs(:,:,ns) = freqAlfaADDnoAvg{ns}';
      allSubjTot(:,:,ns) = [picoAlfaADDnoAvg{ns} freqAlfaADDnoAvg{ns}]';
      group(ns) = 'A';
    else
      for e = 1:size(picoAlfaNOLDnoAvg{ns-16},1) % epochs
        allSubjLabels(e,ns) = 'N';  % all epochs labels
      end
      allSubjPicos(:,:,ns) = picoAlfaNOLDnoAvg{ns-16}';
      allSubjFreqs(:,:,ns) = freqAlfaNOLDnoAvg{ns-16}';
      allSubjTot(:,:,ns) = [picoAlfaNOLDnoAvg{ns-16} freqAlfaNOLDnoAvg{ns-16}]';
      group(ns) = 'N';
    end
    
end

fid = fopen([PATHIN filesep 'ADD (16) x NOLD (19) noAvg classification (10-fold) results.txt'],'wt');
fprintf(fid,'Run\tAccPico\tAccFreq\tAccTot\n');

for r = 1:10
    
 temp = allSubjPicos;
 temp = temp(:,:);
 features_train_tot = temp'; % features in columns
%  size_features_train_tot = size(features_train_tot)
 
 labels_tr = allSubjLabels;
 labels_tr = labels_tr(:);
%  size_labels_tr = size(labels_tr)
%  pause;
 
 grp = unique(labels_tr);
 
 if length(grp) ~= 2
   error('length(grp) ~= 2');
 end
 
 modSVM = fitcsvm(features_train_tot,labels_tr,...
        'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
    
 CV = crossval(modSVM); % 10-fold cross-validation
 
 erroP = kfoldLoss(CV);
 
 temp = allSubjFreqs;
 temp = temp(:,:);
 features_train_tot = temp'; % features in columns
%  size_features_train_tot = size(features_train_tot)
 
 modSVM = fitcsvm(features_train_tot,labels_tr,...
        'Standardize',1,'KernelFunction','rbf','KernelScale','auto');
 
 CV = crossval(modSVM); % 10-fold cross-validation
 
 erroF = kfoldLoss(CV);    
 
 temp = allSubjTot;
 temp = temp(:,:);
 features_train_tot = temp'; % features in columns
%  size_features_train_tot = size(features_train_tot)
 
 modSVM = fitcsvm(features_train_tot,labels_tr,...
        'Standardize',1,'KernelFunction','rbf','KernelScale','auto');
 
 CV = crossval(modSVM); % 10-fold cross-validation
 
 erroT = kfoldLoss(CV);     
 
%fprintf(fid,'Run\tAccPico\tAccFreq\tAccTot\n');
 fprintf(fid,'%d\t%f\t%f\t%f\n',r,100*(1-erroP),100*(1-erroF),100*(1-erroT)); 
 
end
fclose(fid);

