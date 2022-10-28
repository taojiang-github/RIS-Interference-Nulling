close all
% load main_phase_transition_Rician.mat
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

l=plot(num_irs_elements_set,prob_set(:,1),'o-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,2),'d-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,3),'s-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,4),'>-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,5),'v-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,6),'<-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,7),'+-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,8),'*-');hold on
l.MarkerFaceColor = l.Color;
l=plot(num_irs_elements_set,prob_set(:,9),'.-');hold on
l.MarkerFaceColor = l.Color;
grid on
legend('$K=2$','$K=3$','$K=4$','$K=5$','$K=6$','$K=7$','$K=8$','$K=9$','$K=10$','Interpreter','latex','Location','southeast')
xlabel('Number of RIS elements $N$ ','Interpreter','latex')
ylabel('Empirical feasible probability','Interpreter','latex')
xticks(0:20:200)