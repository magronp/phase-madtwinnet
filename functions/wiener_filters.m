% A general function for Wiener filters, under its classical (phase-unaware),
% anisotropic, consistent or consistent anisotropic variants.
%
% Ref: "Consistent anisotropic Wiener filtering for audio source separation",
% Paul Magron, Jonathan Le Roux and Tuomas Virtanen
% IEEE Workshop on Applications of Signal Processing to Audio and Acoustics
% (WASPAA), October 2017
% 
% Inputs:
%     X : F*T STFT of the mixture
%     variances : F*T*J estimates of the sources' variances
%     kappa >= 0 : anisotropy parameter
%     delta >= 0 : consistency parameter
%     Nw : STFT window length
%     hop : hop size of the STFT
%     wien : initial sources estimates (classical or anisotropic Wiener
%     filters)
%
% Outputs:
%     Se : estimated STFT sources

function Se = wiener_filters(X,variances,kappa,delta,Nw,hop,wien)

% Check if a wiener filter is already provided as input for initialization
if nargin<7
    if kappa == 0 % Wiener filter
        wien = variances ./ (sum(variances,3)+eps) .* X;
    else   %anisotropic Wiener filter
        wien = aw_dnn(X,variances,kappa,hop);
    end
end

Se = wien;

% Consistent filter if delta neq 0
if not(delta ==0)
    
    [F,T]=size(X);
    Nfft = 2*(F-1);
    
    voice_wiener = wien(:,:,1);
    mu = angle(wien);
    
    % Anisotropy parameters
    lambda = besseli(1,kappa) ./ besseli(0,kappa);
    rho = (besseli(2,kappa).*besseli(0,kappa) - besseli(1,kappa).^2 )./ besseli(0,kappa).^2;

    % Sources covariance
    gamma = (1-lambda.^2).* variances; gamma_voice = gamma(:,:,1);
    c = rho.*variances .*exp(2*1i*mu);  c_voice = c(:,:,1);

    % Mixture covariance
    gamma_X = sum(gamma,3); c_X = sum(c,3); detGX = gamma_X.^2 - abs(c_X).^2+eps;

    % Posterior covariance
    gamma_post = abs(gamma_voice - ( gamma_X .* (gamma_voice.^2+abs(c_voice).^2) - 2 * gamma_voice .* real(c_voice.*conj(c_X))  )  ./ detGX);
    c_post = c_voice - (2*gamma_voice.*gamma_X.*c_voice - gamma_voice.^2 .* c_X - c_voice.^2 .* conj(c_X)  )  ./ detGX;
    detG = gamma_post.^2-abs(c_post).^2;


    %%% Preconditioned Conjugate Gradient algorithm %%%
    wei=repmat([1; 2*ones(F-2,1); 1],[1 T]);
    se=iSTFT(voice_wiener,Nfft,hop,Nw);
    FSE=voice_wiener-STFT(se,Nfft,hop,Nw);
    r=-delta*FSE;
    z= detG ./ (1+2*delta*(Nw-hop)/Nw*gamma_post + (delta*(Nw-hop)/Nw)^2 * detG) .* ( (gamma_post./(detG+eps) + delta*(Nw-hop)/Nw ).* r + c_post./(detG+eps) .* conj(r) );
    P=z;
    rsold=real(sum(sum(wei.*conj(r).*z)));
    converged=false;

    while ~converged
        p=iSTFT(P,Nfft,hop,Nw);
        FP=P-STFT(p,Nfft,hop,Nw);
        
        AP=( gamma_post .* P - c_post .* conj(P)  )./(detG+eps) +delta*FP;
        alpha=rsold/real(sum(sum(wei.*conj(P).*AP))+realmin);
        voice_wiener=voice_wiener+alpha*P;

        converged=(sum(sum(alpha^2*real(P.*conj(P)))) < 1e-6*sum(sum(real(voice_wiener.*conj(voice_wiener)))));
        r=r-alpha*AP;
        z =  detG ./ (1+2*delta*(Nw-hop)/Nw*gamma_post + (delta*(Nw-hop)/Nw)^2 * detG) .* ( (gamma_post./(detG+eps) + delta*(Nw-hop)/Nw ).* r + c_post./(detG+eps) .* conj(r) );

        rsnew=real(sum(sum(wei.*conj(r).*z)));
        beta=rsnew/(rsold+realmin);
        P=z+beta*P;
        rsold=rsnew;
    end
    
    Se(:,:,1) = voice_wiener;
    Se(:,:,2) = X - voice_wiener;

end

end


% Anisotropic Wiener filter
function Se = aw_dnn(X,variances,kappa,hop)
% Inputs:
%     X : F*T mixture
%     variances : F*T*J variances
%     kappa : concentration parameter
%     hop : STFT overlap (in samples)
%
% Outputs:
%     Se : estimated sources
%     mu : unwrapped phase (sinusoidal model)

[F,T,J] = size(variances);

% Unwrapping Parameters
Nfft = (F-1)*2;
Ct = 2*pi*hop/Nfft;

% Weight parameters
lambda = besseli(1,kappa) / besseli(0,kappa);
rho = (besseli(2,kappa)*besseli(0,kappa) - besseli(1,kappa)^2 )/ besseli(0,kappa)^2;

% Initialization
mu = repmat(angle(X),[1 1 J]);
Se=sqrt(variances) .* exp(1i * mu);

% Loop over time frames
for t=2:T

    % Initialisation with PU
    for j=1:J
        f_inf = get_frequencies_qifft_frame(abs(Se(:,t,j))+eps);
        mu(:,t,j) = angle(Se(:,t-1,j))+Ct*f_inf;
        Se(:,t,j) = abs(Se(:,t,j)) .* exp(1i * mu(:,t,j));
    end

    % Compute prior, means, variances and covariances
    Xtilde = squeeze(Se(:,t,:));

    m = lambda.*Xtilde;
    gamma = (1-lambda.^2).* abs(Xtilde).^2;
    c = rho.*Xtilde.^2;

    m_X = repmat(sum(m,2),[1 J]);
    gamma_X = repmat(sum(gamma,2),[1 J]);        
    c_X = repmat(sum(c,2),[1 J]);

    Xaux = repmat(X(:,t),[1 J]);

    % Get the MMSE estimator
    m_post = m + ( (gamma.*gamma_X - c.*conj(c_X)).*(Xaux-m_X) + (c.*gamma_X - gamma.*c_X).* conj(Xaux-m_X)    ) ./ (gamma_X.^2 - abs(c_X).^2+eps);

    Se(:,t,:) = m_post;

end

end


% Inst. frequency and regions of influence
function [f_inf,f_centr,f_harm] = get_frequencies_qifft_frame(mag)

% Input spectrum expected as row
mag = mag(:)';

%Central peaks
[~,f_centr] = findpeaks(max(mag,10^-6),'MINPEAKHEIGHT',max(mag)*10^-2);

Nfreq = length(f_centr);
f_harm = zeros(1,Nfreq);

if (Nfreq >0)
    % Quadratic interpolation of frequencies
    for ind = 1:Nfreq
        f = f_centr(ind);
        f_harm(ind) = qint(log10(mag(f-1)),log10(mag(f)),log10(mag(f+1)))+f;
    end

    % Frequencies in Regions of influence
    f_inf = zeros(length(mag),1);
    deb = 1;
    index_lim = zeros(1,Nfreq-1);
    
    for ind = 1:(Nfreq-1)
        f = f_centr(ind);
        fp = f_centr(ind+1);
        fin = floor((mag(fp)*f+mag(f)*fp)/(mag(fp)+mag(f)));
        f_inf(deb:fin) = f_harm(ind);
        deb = fin+1;
        index_lim(ind) = fin;
    end

    f_inf(deb:end) = f_harm(end)-1;
    
else
    f_inf = (1:length(mag))'-1;
end

end

% Quadratic Interpolated FFT
function [p,b,a] = qint(ym1,y0,yp1)

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
b = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end
