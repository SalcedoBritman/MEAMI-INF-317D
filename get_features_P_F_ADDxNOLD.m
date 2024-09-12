clear all;
close all;
fclose('all');
clc;

PATHIN = input('Entre com o path do STUDY: ','s');

load([PATHIN filesep 'freqIni.mat']); 
%freqIni : cell array (19 chans)
load([PATHIN filesep 'freqEnd.mat']); 
%freqEnd : cell array (19 chans)

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%[STUDY ALLEEG] = pop_loadstudy('filename', 'ADD x NOLD (16 x 19).study', 'filepath', 'C:\Users\RESEARCH\MESTRADO\FROM 2024\Britman Salcedo (2024)\2024-2 machine learning');
[STUDY ALLEEG] = pop_loadstudy('filename', 'ADD x NOLD (16 x 19).study', 'filepath', PATHIN);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];

% STUDY = pop_statparams(STUDY, 'condstats','on','mcorrect','fdr','alpha',0.05);
% STUDY = pop_specparams(STUDY, 'plotconditions','together','freqrange',[1 45] );
% 
% [STUDY, ~, specfreqs, ~, pgroup] = ...
%       std_specplot(STUDY,ALLEEG,'channels',{EEG(1).chanlocs(1).labels}, 'design', 1);
%   
% chan_freqs_signif = zeros(length(EEG(1).chanlocs),length(pgroup{1}));
% 
% % fid = fopen([PATHIN filesep 'ADDxNOLD_signif_freqs.txt'],'wt');
% % fprintf(fid,'Channel\tSignifFreqRange\n');

numSubj = length(EEG);

picoAlfaADDnoAvg = cell(1,16);
freqAlfaADDnoAvg = cell(1,16);
picoAlfaNOLDnoAvg = cell(1,19);
freqAlfaNOLDnoAvg = cell(1,19);

for ns = 1:numSubj
    
    ns
    
    if ns < 17 % ADD
      picoAlfaADDnoAvg{ns} = zeros(33,19); % epochs x features  
      freqAlfaADDnoAvg{ns} = zeros(33,19); % epochs x features
    else
      picoAlfaNOLDnoAvg{ns-16} = zeros(33,19); % epochs x features
      picoAlfaNOLDnoAvg{ns-16} = zeros(33,19); % epochs x features
    end
    
    CURRENTSTUDY = 0;
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, [1:35] ,'retrieve',ns,'study',1);
      
    psd = zeros(4*EEG.srate+1,EEG.nbchan,33);
    
%     if ns == 1
%         
%       meanADDpsd = zeros(4*EEG.srate+1,EEG.nbchan);
%       meanNOLDpsd = zeros(4*EEG.srate+1,EEG.nbchan);
%       
%     end
    
    for ep = 1:33
        
        epoch = squeeze(EEG.data(:,:,ep))';
        epoch = epoch - ones(size(EEG.data,2),1)*mean(epoch); % elimina nivel DC da epoca
        [psd(:,:,ep),f] = pwelch(epoch,EEG.srate,round(0.9*EEG.srate),8*EEG.srate,EEG.srate);
        f = double(f);
        for ch = 1:EEG.nbchan
          [picoAlfa,indAlfa] = max(psd(f >= 4.0 & f <= 12.0,ch,ep));
          freqAlfa = f(indAlfa) + 4.0;
%           plot(f,psd(:,ch,ep))
%           pause;
          if ns < 17 % ADD
            picoAlfaADDnoAvg{ns}(ep,ch) = 10*log10(picoAlfa);
            freqAlfaADDnoAvg{ns}(ep,ch) = freqAlfa;
          else
            picoAlfaNOLDnoAvg{ns-16}(ep,ch) = 10*log10(picoAlfa);
            freqAlfaNOLDnoAvg{ns-16}(ep,ch) = freqAlfa;
          end
        end

    end
     
%     for ch = 1:EEG.nbchan
%      if ns < 17 % ADD
% %       plot(picoAlfaADD{ns}(:,ch));
%       picoAlfaADD{ns}(:,ch) = filtfilt(ones(1,5),1,picoAlfaADD{ns}(:,ch))/25;
% %       hold on;
% %       plot(picoAlfaADD{ns}(:,ch),'r');
% %       grid on;
% %       title([EEG(1).chanlocs(ch).labels ' - picoAlfaADD' int2str(ns)]);
% %       pause;
% %       hold off;
% %       plot(freqAlfaADD{ns}(:,ch));
%       freqAlfaADD{ns}(:,ch) = filtfilt(ones(1,5),1,freqAlfaADD{ns}(:,ch))/25;
% %       hold on;
% %       plot(freqAlfaADD{ns}(:,ch),'r');
% %       grid on;
% %       title([EEG(1).chanlocs(ch).labels ' - freqAlfaADD' int2str(ns)]);
% %       pause;      
%      else
% %       plot(picoAlfaNOLD{ns-16}(:,ch));
%       picoAlfaNOLD{ns-16}(:,ch) = filtfilt(ones(1,5),1,picoAlfaNOLD{ns-16}(:,ch))/25;
% %       hold on;
% %       plot(picoAlfaNOLD{ns-16}(:,ch),'r');
% %       grid on;
% %       title([EEG(1).chanlocs(ch).labels ' - picoAlfaNOLD' int2str(ns)]);
% %       pause;
% %       hold off;
% %       plot(freqAlfaNOLD{ns-16}(:,ch));
%       freqAlfaNOLD{ns-16}(:,ch) = filtfilt(ones(1,5),1,freqAlfaNOLD{ns-16}(:,ch))/25;
% %       hold on;
% %       plot(freqAlfaNOLD{ns-16}(:,ch),'r');
% %       grid on;
% %       title([EEG(1).chanlocs(ch).labels ' - freqAlfaNOLD' int2str(ns)]);
% %       pause; 
%      end
%     end
    
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
    
end
  
save([PATHIN filesep 'picoAlfaADDnoAvg.mat'],'picoAlfaADDnoAvg'); 
save([PATHIN filesep 'freqAlfaADDnoAvg.mat'],'freqAlfaADDnoAvg');
save([PATHIN filesep 'picoAlfaNOLDnoAvg.mat'],'picoAlfaNOLDnoAvg'); 
save([PATHIN filesep 'freqAlfaNOLDnoAvg.mat'],'freqAlfaNOLDnoAvg');  

% for ch = 1:19
% H = figure;
% plot(f,meanADDpsd(:,ch)/16,'b');
% hold on;
% plot(f,meanNOLDpsd(:,ch)/19,'g');
% title(EEG(1).chanlocs(ch).labels);
% xlabel('f (Hz)');
% ylabel('PSD (dB)');
% legend('ADD','NOLD');
% v = axis;
% axis([1 45 v(3) v(4)]);
% pause;
% saveas(H,[PATHIN filesep ['PSD_' EEG(1).chanlocs(ch).labels '.png']]);
% end


