clear;close all
num_samples = 100;
num_user = 8;
num_elements_irs = 144;
irs_Nh = 12;
Rician_factor = 10;
scale_factor = 120;
sigma2_dB = -10;
file_name = sprintf('./Rician_channel/channel(%d, %d, %d, %d, %d, %d).mat',...
    num_elements_irs, num_user, irs_Nh, num_samples, Rician_factor, scale_factor);
load(file_name);
channel_cascaded = zeros(num_samples,num_user, num_elements_irs, num_user);
for kk=1:num_user
    for jj=1:num_user
        channel_cascaded(:,jj,:,kk)= channel_tx(:, :, jj) .* channel_rx(:, :, kk);
    end
end
num_samples = 100;
loss_set_rand_int = nan(21,num_samples);
loss_set_svd_int = nan(21,num_samples);
for ii=1:num_samples
    A = squeeze(channel_cascaded(ii,:,:,:));
    [v,loss] = projected_gradient(A,1);
    loss_set_svd_int(:,ii) = loss;
    [v,loss] = projected_gradient(A,2);
    loss_set_rand_int(:,ii) = loss;
end
% save main_PG_convergence.mat
%%
% plot(loss_set_rand_int,'Color',[0.00,0.45,0.74])
% hold on
% plot(loss_set_svd_int,'Color','b-')

l=plot(mean(loss_set_rand_int,2),'-V','Color',[0.00,0.45,0.74],'LineWidth',1.5);
l.MarkerFaceColor = l.Color;
hold on
l=plot(mean(loss_set_svd_int,2),'-o','Color',[0.64,0.08,0.18],'LineWidth',1.5);
l.MarkerFaceColor = l.Color;
legend('Random Initialization','Proposed Initialization')
xlabel('Number of iterations','Interpreter','latex')
ylabel('Signal Power $\sum_{k=1}^K |{\bf{a}}_{k,k}^\top \bf{v}|^2$','Interpreter','latex')
grid on
% if ~exist('results', 'dir')
%     mkdir('results');
% end
% file_name = sprintf('./results/main_AP_convergence(%d, %d, %d, %d, %d, %d).mat',...
%     num_elements_irs, num_user, irs_Nh, num_samples, Rician_factor, scale_factor);
% save(file_name,'sigma2_dB', 'num_elements_irs','irs_Nh','num_samples','num_user','Rician_factor','scale_factor');


