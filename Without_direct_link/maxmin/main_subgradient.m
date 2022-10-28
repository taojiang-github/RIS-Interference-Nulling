clear;close all
% rng(1)
%% generate Gaussian channel
num_samples = 500;
num_user_set = [4];
num_irs_elements_set = [16,36];
irs_Nh_set = [4,6];
sigma2_dB_set = 0:-2:-16;
Rician_factor = 10;
scale_factor = 120;
initialization = 0;

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
                if initialization == 0
                    v =exp(1j.*2*pi*rand(num_irs_elements,1));
                elseif initialization == 1
                    v = alter_proj(A);
                elseif initialization == 2
                    v = alter_proj_2(A);
                end
                [v,obj_tmp] = subgradient_maxmin(A,sigma2_dB,v);
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
end
running_time = toc(ts)/60;
fprintf('running time: %3.3f mins\n',toc(ts)/60);
%% save
save main_subgradient1.mat


