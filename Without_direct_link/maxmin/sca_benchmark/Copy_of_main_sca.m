clear;close all
% rng(1)
%% generate Gaussian channel
num_samples = 300;
num_user_set = [4];
num_irs_elements_set = [16];
irs_Nh_set = [4];
sigma2_dB_set = 0:-2:-16;
Rician_factor = 10;
scale_factor = 120;

ts = tic();
rate_set = zeros(length(sigma2_dB_set),length(num_irs_elements_set),length(num_user_set));
rate_set_rand = zeros(length(sigma2_dB_set),length(num_irs_elements_set),length(num_user_set));
for iter_snr=1:length(sigma2_dB_set)
    sigma2_dB = sigma2_dB_set(iter_snr)
    for iter_m = 1:length(num_irs_elements_set)
        num_irs_elements = num_irs_elements_set(iter_m)
        irs_Nh = irs_Nh_set(iter_m);
        for iter_k = 1: length(num_user_set)
            num_user = num_user_set(iter_k)
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
            rate_tmp = 0;
            rate_rnd_tmp = 0;
            parfor ii=1:num_samples
                A = squeeze(channel_cascaded(ii,:,:,:));
                p = ones(num_user,1);
                alpha=ones(num_user,1);
                v =exp(1j.*2*pi*rand(num_irs_elements,1));
                [v, obj_tmp] = sca_func(A,sigma2_dB,v);
                %min rate
                [tmp_rate,tmp_sinr] = compute_minrate(A, v, p, alpha, sigma2_dB);
                rate_tmp = rate_tmp+tmp_rate;
                
                v =exp(1j.*2*pi*rand(num_irs_elements,1));
                v = v./abs(v);
                [tmp_rate_rnd,tmp_sinr_rnd] = compute_minrate(A, v, p, alpha, sigma2_dB);
                rate_rnd_tmp = rate_rnd_tmp+tmp_rate_rnd;
            end
            rate_set(iter_snr,iter_m,iter_k) = rate_tmp/num_samples;
            rate_set_rand(iter_snr,iter_m,iter_k) = rate_rnd_tmp/num_samples;
        end
    end
    save main_sca_16_300.mat
end
running_time = toc(ts)/60;
fprintf('running time: %3.3f mins\n',toc(ts)/60);
%% save
save main_sca_16_300.mat


