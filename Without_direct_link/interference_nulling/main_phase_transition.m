clear;close all
% rng(1)
num_samples = 500;
num_irs_elements_set = 1:200;
num_user_set = 2:10;
irs_Nh_set = ones(length(num_irs_elements_set),1);
Rician_factor = 10;
scale_factor = 120;
sigma2_dB = -10;
channel_model = 'Rician';
% channel_model = 'Sparse';
% channel_model = 'Rayleigh';
num_path = 5;

ts = tic();
signal_power_set = cell(length(num_user_set),length(num_irs_elements_set),num_samples);
interfence_power_set = cell(length(num_user_set),length(num_irs_elements_set),num_samples);
sir_set = cell(length(num_user_set),length(num_irs_elements_set),num_samples);
min_sir_set = nan(length(num_user_set),length(num_irs_elements_set),num_samples);
sum_sir_set = nan(length(num_user_set),length(num_irs_elements_set),num_samples);

for iter_k = 1:length(num_user_set)
    num_user = num_user_set(iter_k)
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
        elseif strcmp(channel_model,'Rayleigh')
            channel_tx = (randn(num_samples,num_irs_elements,num_user)+1j*randn(num_samples,num_irs_elements,num_user))/sqrt(0.5);
            channel_rx = (randn(num_samples,num_irs_elements,num_user)+1j*randn(num_samples,num_irs_elements,num_user))/sqrt(0.5);
        end
        channel_cascaded = zeros(num_samples,num_user, num_irs_elements, num_user);
        for kk=1:num_user
            for jj=1:num_user
                channel_cascaded(:,jj,:,kk)= channel_tx(:, :, jj) .* channel_rx(:, :, kk);
            end
        end
        parfor ii=1:num_samples
            A = squeeze(channel_cascaded(ii,:,:,:));
            v = alter_proj(A);
            [signal_power, interfence_power, sir] = compute_SIR(A, v);
            signal_power_set{iter_k,iter_m,ii} = signal_power;
            interfence_power_set{iter_k,iter_m,ii}  = interfence_power;
            sir_set{iter_k,iter_m,ii}  = sir;
            min_sir_set(iter_k,iter_m,ii) = min(sir);
            sum_sir_set(iter_k,iter_m,ii) = sum(sir);
        end
    end
end
fprintf('running time: %f mins\n',(toc(ts))/60);
if strcmp(channel_model,'Sparse')
    save main_phase_transition_Sparse.mat
elseif strcmp(channel_model,'Rician')
    save main_phase_transition_Rician.mat
elseif strcmp(channel_model,'Rayleigh')
    save main_phase_transition_Rayleigh.mat
end
%%
close all
th = 10^6;
prob_set = zeros(length(num_irs_elements_set),length(num_user_set));
for ii = 1:length(num_irs_elements_set)
    for jj=1:length(num_user_set)
        for kk=1:num_samples
           if min_sir_set(jj,ii,kk)>th
               prob_set(ii,jj) = prob_set(ii,jj)+1;
           end
        end
    end
end
prob_set = prob_set/num_samples;
colormap('gray');   % set colormap
imagesc(flipud(prob_set),[0,1]); 
set(gca, 'ydir', 'reverse')

xlabel('Number of tranceiver pairs $K$ ','Interpreter','latex')
ylabel('Number of IRS elements $N$ ','Interpreter','latex')
yticks(1:10:200)
yticklabels(200:-10:1)
xticks(1:9)
xticklabels(2:10)
grid on
hold on
x = 1:9;
y = 2*x.*(x+1);
plot(x,-y+200,'o-y','MarkerFaceColor','y')
legend('$N=2K(K-1)$','Interpreter','latex','Location','southeast')
colorbar
% axis square 

%%
close all
lb_set = zeros(length(num_user_set),1);
ub_set = zeros(length(num_user_set),1);
for ii = 1:length(num_user_set)
    tmp = find(prob_set(:,ii)==0);
    lb_set(ii) = tmp(end)+1;
    tmp = find(prob_set(:,ii)<=0.95);
    ub_set(ii) = tmp(end)+1;
end
plot(num_user_set,lb_set,'o-')
hold on
plot(num_user_set,ub_set,'s-')
hold on
x = 2:10;
y = 2*x.*(x-1);
plot(x,y,'o-y','MarkerFaceColor','y')
