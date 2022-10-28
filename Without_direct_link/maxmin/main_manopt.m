clear;close all
% rng(1)
%% generate Gaussian channel
num_sample = 100;
num_user = 10;
num_elements_irs_set = [128,256];
sigma2_dB_set = [0,-10];
rate_set = zeros(length(sigma2_dB_set),length(num_elements_irs_set));
rate_set_rand = zeros(length(sigma2_dB_set),length(num_elements_irs_set));

ts = tic();
for iter_snr=1:length(sigma2_dB_set)
    sigma2_dB = sigma2_dB_set(iter_snr);
    for iter_m = 1:length(num_elements_irs_set)
        num_elements_irs = num_elements_irs_set(iter_m);
        channel_tx = 1/sqrt(2)*randn(num_sample,num_elements_irs, num_user)...
            +1/sqrt(2)*1j*randn(num_sample,num_elements_irs, num_user);
        channel_rx = 1/sqrt(2)*randn(num_sample,num_elements_irs, num_user)...
            +1/sqrt(2)*1j*randn(num_sample,num_elements_irs, num_user);
        channel_cascaded = zeros(num_sample,num_user, num_elements_irs, num_user);
        for kk=1:num_user
            for jj=1:num_user
                channel_cascaded(:,jj,:,kk)= channel_tx(:, :, jj) .* channel_rx(:, :, kk);
            end
        end
        rate_tmp = 0;
        rate_rnd_tmp = 0;
        parfor ii=1:num_sample
            A = squeeze(channel_cascaded(ii,:,:,:));
            p = ones(num_user,1);
            alpha=ones(num_user,1);
            v = mamnopt_func(A,p,alpha,sigma2_dB);
            %min rate
            [tmp_rate,tmp_sinr] = compute_minrate(A, v, p, alpha, sigma2_dB);
            rate_tmp = rate_tmp+tmp_rate;
            
             v =exp(1j.*2*pi*rand(num_elements_irs,1));
             v = v./abs(v);
             [tmp_rate_rnd,tmp_sinr_rnd] = compute_minrate(A, v, p, alpha, sigma2_dB);
             rate_rnd_tmp = rate_rnd_tmp+tmp_rate_rnd;
            %sum rate
%             [tmp_rate,tmp_interference] = compute_rate(A, v, p, alpha, sigma2_dB);
%             interference_set(ii,iter_m) = tmp_interference;
%             rate_tmp = rate_tmp+tmp_rate;
%             prob_tmp = prob_tmp+(tmp_interference<1e-4)*1.0;
            
        end
        rate_set(iter_snr,iter_m) = rate_tmp/num_sample;
        rate_set_rand(iter_snr,iter_m) = rate_rnd_tmp/num_sample;
    end    
end
running_time = toc(ts)/60;
fprintf('running time: %3.3f mins\n',toc(ts)/60);

if ~exist('results', 'dir')
    mkdir('results');
end
file_name = sprintf('./results/Manopt_Gaussian_k%d_%d.mat',num_user,ts)
save(file_name,'sigma2_dB_set', 'num_elements_irs_set', 'rate_set','rate_set_rand','running_time');


