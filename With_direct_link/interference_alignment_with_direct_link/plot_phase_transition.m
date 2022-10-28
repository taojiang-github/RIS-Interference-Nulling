close all
load main_phase_transition_Rician.mat
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

lb_set = zeros(length(num_user_set),1);
ub_set_95 = zeros(length(num_user_set),1);
for ii = 1:length(num_user_set)
    tmp = find(prob_set(:,ii)==0);
    lb_set(ii) = tmp(end)+1;
    tmp = find(prob_set(:,ii)<=0.95);
    ub_set_95(ii) = tmp(end)+1;
end
x = 2:10;
y = 2*x.*(x-1);
plot(x,y,'-*')
hold on
plot(num_user_set,lb_set,'o-')
hold on
plot(num_user_set,ub_set_95,'s-')
hold on

load main_phase_transition_Sparse.mat
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

lb_set = zeros(length(num_user_set),1);
ub_set_95 = zeros(length(num_user_set),1);
for ii = 1:length(num_user_set)
    tmp = find(prob_set(:,ii)==0);
    lb_set(ii) = tmp(end)+1;
    tmp = find(prob_set(:,ii)<=0.95);
    ub_set_95(ii) = tmp(end)+1;
end
plot(num_user_set,lb_set,'o--')
hold on
plot(num_user_set,ub_set_95,'s--')
hold on
legend('$N=2K(K-1)$','0\% success (Rician)','95\% success (Rician)','0\% success (Sparse)','95\% success (Sparse)','Interpreter','latex','Location','southeast')
xlabel('Number of tranceiver pairs $K$ ','Interpreter','latex')
ylabel('Number of RIS elements $N$ ','Interpreter','latex')
grid on
yticks(0:10:200)