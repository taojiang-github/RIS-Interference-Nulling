clear;close all
ts = tic();
%% paramters
num_samples = 500;
num_user = 6;
num_irs_elements_set = 40:10:100;
irs_Nh = 10;
sigma2_dB = -15;
Rician_factor = 10;
scale_factor = 120;

rate_set1 = zeros(length(num_irs_elements_set),1);
rate_set2 = zeros(length(num_irs_elements_set),1);
rate_set3 = zeros(length(num_irs_elements_set),1);
rate_set4 = zeros(length(num_irs_elements_set),1);
rate_set5 = zeros(length(num_irs_elements_set),1);
rate_set_rand = zeros(length(num_irs_elements_set),1);

for iter_m=1:length(num_irs_elements_set)
    num_irs_elements = num_irs_elements_set(iter_m);
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
    rate_tmp1 = 0; rate_tmp2 = 0; rate_tmp3 = 0;rate_tmp4 = 0; rate_tmp5 = 0;
    rate_rnd_tmp = 0;
    parfor ii=1:num_samples
        A = squeeze(channel_cascaded(ii,:,:,:));
        p = ones(num_user,1);
        alpha=ones(num_user,1);
        
        % random initialization
        v0 =exp(1j.*2*pi*rand(num_irs_elements,1));
        v = mamnopt_func(A,p,alpha,sigma2_dB,v0);
        rate_tmp1 = rate_tmp1+compute_rate(A, v, p, alpha, sigma2_dB);
        
        % zf initialization
        v0 = alter_proj(A);
        rate_tmp4 = rate_tmp4+compute_rate(A, v0, p, alpha, sigma2_dB);
        v = mamnopt_func(A,p,alpha,sigma2_dB,v0);
        rate_tmp2 = rate_tmp2+compute_rate(A, v, p, alpha, sigma2_dB);
        
        % zf initialization
        v0 = alter_proj_2(A);
        rate_tmp5 = rate_tmp5+compute_rate(A, v0, p, alpha, sigma2_dB);
        v = mamnopt_func(A,p,alpha,sigma2_dB,v0);
        rate_tmp3 = rate_tmp3+compute_rate(A, v, p, alpha, sigma2_dB);
        
        % random irs
        v =exp(1j.*2*pi*rand(num_irs_elements,1));
        v = v./abs(v);
        rate_rnd_tmp = rate_rnd_tmp+compute_rate(A, v, p, alpha, sigma2_dB);
    end
    rate_set1(iter_m) = rate_tmp1/num_samples;
    rate_set2(iter_m) = rate_tmp2/num_samples;
    rate_set3(iter_m) = rate_tmp3/num_samples;
    rate_set4(iter_m) = rate_tmp4/num_samples;
    rate_set5(iter_m) = rate_tmp5/num_samples;
    rate_set_rand(iter_m) = rate_rnd_tmp/num_samples;
end
running_time = toc(ts)/60;
fprintf('running time: %3.3f mins\n',toc(ts)/60);
save main_rate_N_manopt2.mat
%% plot
% plot(num_irs_elements_set, rate_set1, 'o-'); hold on
% plot(num_irs_elements_set, rate_set2, 's-'); hold on
% plot(num_irs_elements_set, rate_set3, 'v-'); hold on
% plot(num_irs_elements_set, rate_set4, '<-'); hold on
% plot(num_irs_elements_set, rate_set5, '>-'); hold on
% legend('Random + Manopt','Random + ZF + Manopt','Eig + ZF + Manopt', 'ZF', 'Eig + ZF')
% xlabel('Number of IRS elements $N$ ','Interpreter','latex')
% ylabel('Sum rate (bps/Hz) ','Interpreter','latex')

l=plot(num_irs_elements_set, rate_set4, '<-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set, rate_set5, '>-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set, rate_set1, 'o-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set, rate_set2, 's-'); hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set, rate_set3, 'v-'); hold on
l.MarkerFaceColor = l.Color;
legend('Random Init. + AP','Eig. Init. + AP', 'Random Init. + RCG',... 
'Random Init. + AP + RCG', 'Eig. Init. + AP + RCG','Location','southeast','Interpreter','latex')
xlabel('Number of RIS elements $N$ ','Interpreter','latex')
ylabel('Sum rate (bps/Hz) ','Interpreter','latex')
grid on
