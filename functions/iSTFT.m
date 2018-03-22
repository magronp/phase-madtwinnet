function x = iSTFT(X, Nfft, hop, Nw, wtype)

if nargin<5
    wtype = 'hamming';
end

[F,T,J] = size(X);

% Synthesis window
win = gl_wind(Nw,hop,Nfft,wtype);

% Time-domain signal
hw_1 = floor((Nw+1)/2);
hw_2 = floor(Nw/2);
xlength = (T-1)*hop + Nw;
x = zeros(xlength,J);

% Overalp add
n_in = 1;
for t = 1:T
    
    % Reshape the spectrum so it's Hermitian
    spectro = zeros(Nfft,J);
    spectro(1:F,:) = squeeze(X(:,t,:));
    spectro(F+1:end,:) = conj(squeeze(X(end-1:-1:2,t,:)));
    ifft_buffer = ifft(spectro);
    
    % Roll-back the 0-phase windowing
    xs = zeros(Nw,J);
    xs(1:hw_2,:) = ifft_buffer(end-hw_2+1:end,:);
    xs(hw_2+1:end,:) = ifft_buffer(1:hw_1,:);
    
    % Overlap-add
    time = n_in:(n_in+Nw-1);
    x(time,:) = x(time,:) + win .* xs;
    
    n_in = n_in+hop;
end

% Remove extra zeros added at analysis
x = x((3*hop+1):(xlength-3*hop),:)';


end


function syn_w = gl_wind(Nw,hop,Nfft,wtype)

switch wtype
    case 'hann'
        syn_w = hann(Nw)/sqrt(Nfft);
    case 'hamming'
        syn_w = hamming(Nw)/sqrt(Nfft);
    case 'blackman'
        syn_w = blackman(Nw)/sqrt(Nfft);
end

syn_w_prod = syn_w.^2;
redundancy = floor(Nw/hop);
env = zeros(Nw,1);

for k=-redundancy:redundancy+1

    env_ind = hop * k;
    win_ind = 1:Nw;
    env_ind = env_ind + win_ind;

    valid = ((env_ind > 0) & (env_ind <= Nw));

    env_ind = env_ind(valid) ;
    win_ind = win_ind(valid) ;
    env(env_ind) = env(env_ind) + syn_w_prod(win_ind);

end

syn_w = syn_w./env;
  
end
