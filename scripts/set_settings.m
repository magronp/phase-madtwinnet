% Set the settings used in the experiments

Fs = 44100;
Nsongs = 50;

% STFT parameters
Nfft = 4096;
Nw = 2049;
hop = 384;

context_length = 10;

% Paths
dataset_path = 'datasets/';
magnitudes_path = 'magnitude_spectrograms/';
audio_path = 'audio_files/';
metrics_path = 'metrics/';

% Algorithms
algos = {'Baseline','PU-Iter','CAW'};
Nalgo = length(algos);
Npuiter = 50;
kappa_caw = 0.01;
delta_caw = 0.1;