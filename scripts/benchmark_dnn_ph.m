clear all; close all; clc;
set_settings;

% Separation
for ind=1:Nsongs
    
    clc; fprintf('song %d / %d \n',ind,Nsongs);
    
    %%% Get the data
    [sm,Sm,X,V_estim] = getdata_singing_sep(dataset_path,magnitudes_path,ind,Nfft,hop,Nw,context_length);
    [F,T,J] = size(Sm);
    variances = V_estim.^2;
    
    
    %%% Phase retrieval
    S_estim = zeros(F,T,J,Nalgo);

    % Baseline
    bl = V_estim .* exp(1i*angle(X));
    S_estim(:,:,:,1) = bl;
    
    % PU-Iter
    S_estim(:,:,:,2) = pu_iter(X,bl,hop,Npuiter);  
    
    % Consistent Anisotropic Wiener filter
    S_estim(:,:,:,3) = wiener_filters(X,variances,kappa_caw,delta_caw,Nw,hop);
    
    
    %%% Synthesis and record
    
    % Time-domain synthesis
    s_estim = zeros(J,length(sm),Nalgo);
    for al=1:Nalgo
        s_estim(:,:,al) = real(iSTFT(S_estim(:,:,:,al),Nfft,hop,Nw));
    end
    
    % Remove samples for which the estimation is irelevant (MaD TwinNet)
    s_estim = s_estim(:,context_length*hop+1:end,:);
    sm = sm(:,context_length*hop+1:end);
    
    % Record
    audiowrite(strcat(audio_path,int2str(ind),'_voice_orig.wav'),sm(1,:),Fs);
    audiowrite(strcat(audio_path,int2str(ind),'_music_orig.wav'),sm(2,:),Fs);
    for al = 1:Nalgo
        audiowrite(strcat(audio_path,int2str(ind),'_voice_estim_',int2str(al),'.wav'),squeeze(s_estim(1,:,al)),Fs);
        audiowrite(strcat(audio_path,int2str(ind),'_music_estim_',int2str(al),'.wav'),squeeze(s_estim(2,:,al)),Fs);
    end
    
end


% Compute BSS score
SDR = []; SIR = []; SAR = [];

for ind=1:Nsongs
    
    clc; fprintf('BSS - Song %d / %d \n',ind,Nsongs);
     
    % Original Files path
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
sdaux = SDR; li = isnan(sdaux(1,:)); sdaux(:,li) = [];
siaux = SIR; li = isnan(siaux(1,:)); siaux(:,li) = [];
saaux = SAR; li = isnan(saaux(1,:)); saaux(:,li) = [];

save(strcat(metrics_path,'bss_phase-twinnet.mat'),'SDR','SIR','SAR');

