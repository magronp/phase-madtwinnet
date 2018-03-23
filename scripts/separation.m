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
