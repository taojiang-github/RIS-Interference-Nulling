clear;close all
% rng(1)
num_samples = 2000;
path_loss_direct_set = [0,0.2,0.4,0.6,0.8,1];
num_irs_elements_set = 80:130;
num_user = 7;
irs_Nh_set = ones(length(num_irs_elements_set),1);
Rician_factor = 10;
scale_factor = 120;
sigma2_dB = -10;
channel_model = 'Rician';
% channel_model = 'Sparse';
% num_path = 5;

ts = tic();
signal_power_set = cell(length(path_loss_direct_set),length(num_irs_elements_set),num_samples);
interfence_power_set = cell(length(path_loss_direct_set),length(num_irs_elements_set),num_samples);
sir_set = cell(length(path_loss_direct_set),length(num_irs_elements_set),num_samples);
min_sir_set = nan(length(path_loss_direct_set),length(num_irs_elements_set),num_samples);
sum_sir_set = nan(length(path_loss_direct_set),length(num_irs_elements_set),num_samples);
eta_set = nan(length(path_loss_direct_set),length(num_irs_elements_set));
for iter_k = 1:length(path_loss_direct_set)
    path_loss_direct = path_loss_direct_set(iter_k)
    for iter_m = 1:length(num_irs_elements_set)
        num_irs_elements = num_irs_elements_set(iter_m);
        if num_irs_elements<num_user*(num_user-1)
           continue; 
        end
        irs_Nh = irs_Nh_set(iter_m);
        if strcmp(channel_model,'Sparse')
            [channel_tx,~]=generate_channel_Sparse(num_irs_elements, num_user, irs_Nh, ...
                num_samples,num_path);
            [channel_rx,~]=generate_channel_Sparse(num_irs_elements, num_user, irs_Nh, ...
                num_samples,num_path);
        elseif strcmp(channel_model,'Rician')
            [channel_tx,~]=generate_channel_Rician(num_irs_elements, num_user, irs_Nh, ...
                num_samples, Rician_factor, scale_factor,1);
            [channel_rx,~]=generate_channel_Rician(num_irs_elements, num_user, irs_Nh, ...
                num_samples, Rician_factor, scale_factor,2);
        end
        channel_cascaded = zeros(num_samples,num_user, num_irs_elements, num_user);
        for kk=1:num_user
            for jj=1:num_user
                channel_cascaded(:,jj,:,kk)= channel_tx(:, :, jj) .* channel_rx(:, :, kk);
            end
        end
        channel_direct = randn(num_samples,num_user,num_user)+1j*randn(num_samples,num_user,num_user);
        channel_direct = path_loss_direct * sqrt(0.5)*channel_direct;
        eta = 0;
        parfor ii=1:num_samples
            A = squeeze(channel_cascaded(ii,:,:,:));
            B = squeeze(channel_direct(ii,:,:));
            eta = eta+compute_eta(A,B);
            v = alter_proj(A,B);
            [signal_power, interfence_power, sir] = compute_SIR(A,B, v);
            signal_power_set{iter_k,iter_m,ii} = signal_power;
            interfence_power_set{iter_k,iter_m,ii}  = interfence_power;
            sir_set{iter_k,iter_m,ii}  = sir;
            min_sir_set(iter_k,iter_m,ii) = min(sir);
            sum_sir_set(iter_k,iter_m,ii) = sum(sir);
        end
        eta_set(iter_k,iter_m) = eta/num_samples;
    end
end
fprintf('running time: %f mins\n',(toc(ts))/60);
if strcmp(channel_model,'Sparse')
    save main_phase_transition2_Sparse.mat
elseif strcmp(channel_model,'Rician')
    save main_phase_transition2_Rician.mat
end
%%
close all
th = 10^6;
prob_set = zeros(length(num_irs_elements_set),length(path_loss_direct_set));
for ii = 1:length(num_irs_elements_set)
    for jj=1:length(path_loss_direct_set)
        for kk=1:num_samples
           if min_sir_set(jj,ii,kk)>th
               prob_set(ii,jj) = prob_set(ii,jj)+1;
           end
        end
    end
end
prob_set = prob_set/num_samples;
l=plot(num_irs_elements_set,prob_set(:,1),'o-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,2),'d-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,3),'d-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,4),'d-');hold on
l.MarkerFaceColor = l.Color;
% xlabel('$\eta$ ','Interpreter','latex')
xlabel('Number of RIS elements $N$ ','Interpreter','latex')
ylabel('Empirical feasible probability','Interpreter','latex')
legend('$p=0$','$p=0.1$','$p=0.5$','$p=1$','Interpreter','latex','Location','southeast')
grid on
hold on


