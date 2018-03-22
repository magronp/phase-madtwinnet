%  Load data in Matlab for the phase recovery task
%
% 
% Inputs :
%     dataset_path : path of the dataset
%     magnitudes_path : path of the estimated magnitudes
%     ind : song index (between 0 and 49 for DSD100)
%     Nfft : number of FFT points
%     hop : hop size of the STFT
%     Nw : STFT window length
%     context_length : number of frames of the context (in MaDTwinnet)
% 
% Outputs :
%     sm : time-domain original sources
%     Sm : STFT of the original sources
%     X : STFT of the original mixture
%     V_estim : estimated magnitude spectrograms

function [sm,Sm,X,V_estim] = getdata_singing_sep(dataset_path,magnitudes_path,ind,Nfft,hop,Nw,context_length)

%%% True sources

% Time-domain sources
sources_path = strcat(dataset_path,'Sources/Test/');
aux = dir(sources_path);
song_name = aux(ind+2).name;

bb = audioread(strcat(sources_path,song_name,'/bass.wav'));
dd = audioread(strcat(sources_path,song_name,'/drums.wav'));
oo = audioread(strcat(sources_path,song_name,'/other.wav'));
vv = audioread(strcat(sources_path,song_name,'/vocals.wav'));
mus = mean(bb+dd+oo,2)'; voc = mean(vv,2)';
sm = [voc;mus]+eps;

% STFT and remove initial context frames
Sm = STFT(sm,Nfft,hop,Nw);
Sm = Sm(:,context_length+1:end,:);


%%% Estimated voice magnitude

% song index
aux = int2str(ind-1);
if ind<=10
    index_song = strcat('0',aux);
else
    index_song = aux;
end

% load V
load(strcat(magnitudes_path,index_song,'_source_voice_predicted_spectrogram.mat'));
aux = x; [T1,T2,F]=size(aux);
aux = permute(aux,[2 1 3]);
aux = reshape(aux,[T2*T1 1 F]);
V_voice_estim = squeeze(aux)';



%%% Reshape the STFTs so that orig and estim have the same size
T = min(size(V_voice_estim,2),size(Sm,2));
V_voice_estim = V_voice_estim(:,1:T);
Sm = Sm(:,1:T,:);
X = sum(Sm,3);
Vx = abs(X);
V_estim = zeros(size(Sm)); V_estim(:,:,1) = V_voice_estim; V_estim(:,:,2) = Vx-V_voice_estim;

% Time domain original sources (with the good size)
sm = iSTFT(Sm,Nfft,hop,Nw);

end