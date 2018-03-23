clear all; close all; clc;
set_settings;

SDR = []; SIR = []; SAR = [];

for ind=1:Nsongs
    
    clc; fprintf('BSS - Song %d / %d \n',ind,Nsongs);
     
    % Original Files
    s1 = audioread(strcat(audio_path,int2str(ind),'_voice_orig.wav'));
    s2 = audioread(strcat(audio_path,int2str(ind),'_music_orig.wav'));
    sm = [s1 s2]';
    
    % Loop over the algorithms
    sd_aux = []; si_aux = []; sa_aux = [];
   for al=1:Nalgo
       s_estim1 =  audioread(strcat(audio_path,int2str(ind),'_voice_estim_',int2str(al),'.wav'));
       s_estim2 =  audioread(strcat(audio_path,int2str(ind),'_music_estim_',int2str(al),'.wav'));
       se = [s_estim1 s_estim2]';
       [sdr,~,sir,sar] = bss_eval_images_framewise(se,sm);
       sd_aux = [sd_aux ; sdr(1,:)];
       si_aux = [si_aux ; sir(1,:)];
       sa_aux = [sa_aux ; sar(1,:)];
   end
   
   SDR = [SDR sd_aux]; SIR = [SIR si_aux]; SAR = [SAR sa_aux];
   
end

% Remove NaN
sdaux = SDR; li = isnan(sdaux(1,:));
SDR(:,li) = []; SIR(:,li) = []; SAR(:,li) = [];
save(strcat(metrics_path,'bss_phase-twinnet.mat'),'SDR','SIR','SAR');