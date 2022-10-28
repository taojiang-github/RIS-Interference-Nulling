clear;close all
ts = tic();
%% paramters
num_samples = 500;
num_user = 8;
num_irs_elements = 144;
irs_Nh = 12;
sigma2_dB_set = 0:-2:-16;
% sigma2_dB_set = [-6,-10];
power_set = -sigma2_dB_set;
Rician_factor = 10;
scale_factor = 120;
%% generate channels
[channel_tx,~]=generate_channel_Rician(num_irs_elements, num_user, irs_Nh, ...
    num_samples, Rician_factor, scale_factor,1);
[channel_rx,~]=generate_channel_Rician(num_irs_elements, num_user, irs_Nh, ...
    num_samples, Rician_factor, scale_factor,2);
channel_cascaded = zeros(num_samples,num_user, num_irs_elements, num_user);
for kk=1:num_user
    for jj=1:num_user
        channel_cascaded(:,jj,:,kk)= channel_tx(:, :, jj) .* channel_rx(:, :, kk);
    end
end
%%
rate_set1 = zeros(length(sigma2_dB_set),1);
rate_set2 = zeros(length(sigma2_dB_set),1);
rate_set3 = zeros(length(sigma2_dB_set),1);
rate_set4 = zeros(length(sigma2_dB_set),1);
rate_set5 = zeros(length(sigma2_dB_set),1);
rate_set_rand = zeros(length(sigma2_dB_set),1);
for iter_snr=1:length(sigma2_dB_set)
    sigma2_dB = sigma2_dB_set(iter_snr)
    rate_tmp1 = 0; rate_tmp2 = 0; rate_tmp3 = 0;rate_tmp4 = 0; rate_tmp5 = 0;
    rate_rnd_tmp = 0;
    parfor ii=1:num_samples
        A = squeeze(channel_cascaded(ii,:,:,:));
        p = ones(num_user,1);
        alpha=ones(num_user,1);
        
%         random initialization
        v0 =exp(1j.*2*pi*rand(num_irs_elements,1));
        v = fp_sca_func(A,p,alpha,sigma2_dB,v0);
        rate_tmp1 = rate_tmp1+compute_rate(A, v, p, alpha, sigma2_dB);
        
        % zf initialization
        v0 = alter_proj(A);
        rate_tmp4 = rate_tmp4+compute_rate(A, v0, p, alpha, sigma2_dB);
        v = fp_sca_func(A,p,alpha,sigma2_dB,v0);
        rate_tmp2 = rate_tmp2+compute_rate(A, v, p, alpha, sigma2_dB);
%         
%         % zf initialization
        v0 = alter_proj_2(A);
        rate_tmp5 = rate_tmp5+compute_rate(A, v0, p, alpha, sigma2_dB);
        v = fp_sca_func(A,p,alpha,sigma2_dB,v0);
        rate_tmp3 = rate_tmp3+compute_rate(A, v, p, alpha, sigma2_dB);
%         
        % random irs
        v =exp(1j.*2*pi*rand(num_irs_elements,1));
        v = v./abs(v);
        rate_rnd_tmp = rate_rnd_tmp+compute_rate(A, v, p, alpha, sigma2_dB);
    end
    rate_set1(iter_snr) = rate_tmp1/num_samples;
    rate_set2(iter_snr) = rate_tmp2/num_samples;
    rate_set3(iter_snr) = rate_tmp3/num_samples;
    rate_set4(iter_snr) = rate_tmp4/num_samples;
    rate_set5(iter_snr) = rate_tmp5/num_samples;
    rate_set_rand(iter_snr) = rate_rnd_tmp/num_samples;
    save main_rate_SNR_fp_sca.mat
end
running_time = toc(ts)/60;
fprintf('running time: %3.3f mins\n',toc(ts)/60);
save main_rate_SNR_fp_sca.mat
%% plot
l=plot(power_set, rate_set4, '<-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(power_set, rate_set5, '>-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(power_set, rate_set1, 'o-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(power_set, rate_set2, 's-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(power_set, rate_set3, 'v-'); hold on
l.MarkerFaceColor = l.Color;
xlim([0,16]);
legend('Random Init. + AP','Eig. Init. + AP', 'Random Init. + FP (SCA)',... 
'Random Init. + AP + FP (SCA)', 'Eig. Init. + AP + FP (SCA)','Location','southeast','Interpreter','latex')
xlabel('Transimit power (dBm) ','Interpreter','latex')
ylabel('Sum rate (bps/Hz) ','Interpreter','latex')
grid on
