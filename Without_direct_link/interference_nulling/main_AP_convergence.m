clear;close all
num_samples = 100;
num_user = 8;
num_irs_elements = 144;
irs_Nh = 12;
Rician_factor = 10;
scale_factor = 120;
sigma2_dB = -10;
% file_name = sprintf('./Rician_channel/channel(%d, %d, %d, %d, %d, %d).mat',...
%     num_irs_elements, num_user, irs_Nh, num_samples, Rician_factor, scale_factor);
% load(file_name);

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
loss_set = nan(601,num_samples);
loss_set_pgd = nan(6001,num_samples);
for ii=1:num_samples
    A = squeeze(channel_cascaded(ii,:,:,:));
    [v,loss] = alter_proj(A);
    loss_set(:,ii) = loss.interference2;
    [v,loss] = proj_gradient(A, 1);
    loss_set_pgd(:,ii) = loss.interference2;
end
%%
semilogy(loss_set,'Color',[0.00,0.45,0.74]); hold on
semilogy(loss_set_pgd,'Color','k'); hold on
xlabel('Number of iterations','Interpreter','latex')
% ylabel('Interference $\|\bf{{A}}^\top\bf{{v}}\|_2^2$','Interpreter','latex')
ylabel('Maximum intereference to signal ratio','Interpreter','latex')
grid on
% save main_AP_convergence.mat
% if ~exist('results', 'dir')
%     mkdir('results');
% end
% file_name = sprintf('./results/main_AP_convergence(%d, %d, %d, %d, %d, %d).mat',...
%     num_irs_elements, num_user, irs_Nh, num_samples, Rician_factor, scale_factor);
% save(file_name,'sigma2_dB', 'num_irs_elements','irs_Nh','num_samples','num_user','Rician_factor','scale_factor');


