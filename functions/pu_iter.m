% Iterative source separation procedure using the phase unwrapping technique
% initialization.
%
% Ref:
% "Model-based STFT phase recovery for audio source separation",
% Paul Magron, Roland Badeau and Bertrand David,
% IEEE/ACM Transactions on Speech, Audio, and Language Processing
% Vol 26, number 6, June 2018
%
% Inputs:
%     X : F*T mixture
%     Xe : F*T*J initial sources
%     hop : STFT overlap (in samples)
%     Nit : number of iterations
%
% Outputs:
%     Ye : estimated sources


function [Ye] = pu_iter(X,Xe,hop,Nit)

% Parameters
[F,T,J] = size(Xe);
Nfft = (F-1)*2;
phase_increment = 2*pi*hop/Nfft;

% Weights: Wiener gain
G = abs(Xe).^2./(repmat(sum(abs(Xe).^2,3),[1 1 J])+eps);

% Initial values
Ye=Xe;

% Loop over time frames
for t=2:T
    
        % Initialization
        for j=1:J
            f_inf = get_frequencies_qifft_frame(abs(Ye(:,t,j))+eps);
            phiaux = angle(Ye(:,t-1,j))+phase_increment*f_inf;
            Ye(:,t,j) = abs(Ye(:,t,j)) .* exp(1i * phiaux);
        end
        
        Yaux = squeeze(Ye(:,t,:));
        Gaux = squeeze(G(:,t,:));
            
        % Iiterative procedure
        E = X(:,t) - sum(Yaux,2);
        v = abs(Yaux);
        for iter =1:Nit
            Yaux = Yaux + repmat(E,[1 J]).*Gaux;
            Yaux = Yaux ./ (abs(Yaux)+eps) .* v;
            E = X(:,t) - sum(Yaux,2);
        end
        
        % Update the sources in frame t
        Ye(:,t,:) = Yaux;
        
end

end


% PU: frequencies and regions of influence
function [f_inf,f_centr,f_harm] = get_frequencies_qifft_frame(v)

v = v(:)';

%Central peaks
%[~,f_centr] = findpeaks(v,'MINPEAKHEIGHT',0.01*max(v));
[~,f_centr] = findpeaks(max(v,10^-6),'MINPEAKHEIGHT',max(v)*10^-2);

Nfreq = length(f_centr);
f_harm = zeros(1,Nfreq);

if (Nfreq >0)
    % QIFFT
    for ind = 1:Nfreq
        f = f_centr(ind);
        f_harm(ind) = qint(log10(v(f-1)),log10(v(f)),log10(v(f+1)))+f;
    end

    % Regions of influence
    f_inf = zeros(length(v),1);
    deb = 1;
    index_lim = zeros(1,Nfreq-1);
    
    for ind = 1:(Nfreq-1)
        f = f_centr(ind);
        fp = f_centr(ind+1);
        fin = floor((v(fp)*f+v(f)*fp)/(v(fp)+v(f)));
        f_inf(deb:fin) = f_harm(ind);
        deb = fin+1;
        index_lim(ind) = fin;
    end

    f_inf(deb:end) = f_harm(end);
    
else
    f_inf = (1:length(v))';
end

%  Frequencies start from 0
f_inf = f_inf-1;

end

% Quadratic Interpolated FFT
function [p,b,a] = qint(ym1,y0,yp1)

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
b = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end
