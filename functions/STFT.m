function X = STFT(y, Nfft, hop, Nw, wtype)

if nargin<5
    wtype = 'hamming';
end

% y expected as column vector for each source, so L*J
y = y';
J = size(y,2);
y = [zeros(3*hop,J) ; y ; zeros(3*hop,J)];

% Window
switch wtype
    case 'hann'
        w = hann(Nw)/sqrt(Nfft);
    case 'hamming'
        w = hamming(Nw)/sqrt(Nfft);
    case 'blackman'
        w = blackman(Nw)/sqrt(Nfft);
end
win = repmat(w,[1 J]);


% STFT
%T = floor((size(y,1))/hop);     %kostas style : too much extra info.
T = floor((size(y,1)-Nw)/hop)+1;    %my style: allows to have constant size after 1 pass
X = zeros(Nfft/2+1,T,J);

t=1;
n_in = 1;
n_end = length(y)-Nw+1;

hw_1 = floor((Nw+1)/2);
hw_2 = floor(Nw/2);

while n_in<=n_end
    time = 1 + (t-1) * hop : Nw + (t-1) * hop;
    xw = y(time,:) .* win;
    
    fft_buffer = zeros(Nfft,J);
    fft_buffer(1:hw_1,:) = xw(hw_2+1:end,:);
    fft_buffer(end-hw_2+1:end,:) = xw(1:hw_2,:);
    
    aux = fft(fft_buffer,Nfft);
    X(:,t,:) = aux(1:Nfft/2+1,:);
    t =t+1;
    n_in = n_in+hop;
end
    
end