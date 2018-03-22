clear all; close all; clc;
set_settings;

Kappa = [0 10.^(-3:1)]; Nk = length(Kappa);
Delta = [0 10.^(-3:1)]; Nd = length(Delta);
score = zeros(Nk,3,Nsongs);

for ind=1:Nsongs
    
    clc; fprintf('song %d / %d \n',ind,Nsongs);

    % Get the data: orig mixture and estimated variances
    [sm,Sm,X,V_estim] = getdata_singing_sep(dataset_path,magnitudes_path,ind,Nfft,hop,Nw,context_length);
    sm = sm(:,context_length*hop+1:end);
    [F,T,J] = size(Sm);
    variances = V_estim.^2;
    
    % Initialize the score matrix
    L=length(sm); nwin = floor( (L-30*Fs + 15*Fs) / (15*Fs));
    score =zeros(3,nwin,Nk,Nd);
    
    % Loop over Kappa
    for k=1:Nk
        fprintf('kappa %d / %d \n',k,Nk);
        kappa = Kappa(k);
        anis_wien = wiener_filters(X,variances,kappa,0,Nw,hop);
        
        % Loop over Delta
        for d=1:Nd
            fprintf('delta %d / %d \n',d,Nd);
            delta = Delta(d);
            
            % Consistent Anisotropic Wiener fitlering
            cons_anis_wien = wiener_filters(X,variances,kappa,delta,Nw,hop,anis_wien);
             
            % Time-domain synthesis
            se = real(iSTFT(cons_anis_wien,Nfft,hop,Nw));
            
            % Remove samples for which the estimation is irelevant (MaD TwinNet)
            se = se(:,context_length*hop+1:end);

            % Score
            [sdr,~,sir,sar] = bss_eval_images_framewise(se,sm);
            score(:,:,k,d) = [sdr(1,:);sir(1,:);sar(1,:)];
        end
    end
    
    % Save score for the current song
    save(strcat(metrics_path,'caw/',int2str(ind),'.mat'),'score');
    
end


% Get the overall score

score_total = [];
for ind=[1:35 37:44 46:50]
    load(strcat(metrics_path,int2str(ind),'.mat'));
    score_total = cat(2,score_total,score);
end
li = isnan(score_total(1,:,3,3));
wo = score_total; wo(:,li,:,:) = [];

po = squeeze(median(wo,2));
sdr = squeeze(po(1,1:end-1,:));
sir = squeeze(po(2,1:end-1,:));
sar = squeeze(po(3,1:end-1,:));

figure;
gcmap = repmat(1 - ((0:63)'/64),[1 3]);
colormap(gcmap);
subplot(1,3,1); imagesc(sdr); axis xy; set(gca,'xticklabel',Delta,'yticklabel',Kappa(1:end-1)); xlabel('\delta','fontsize',16); ylabel('\kappa','fontsize',16); title('SDR');
subplot(1,3,2); imagesc(sir); axis xy; set(gca,'xticklabel',Delta,'yticklabel',Kappa(1:end-1)); xlabel('\delta','fontsize',16); title('SIR');
subplot(1,3,3); imagesc(sar); axis xy; set(gca,'xticklabel',Delta,'yticklabel',Kappa(1:end-1)); xlabel('\delta','fontsize',16); title('SAR');
