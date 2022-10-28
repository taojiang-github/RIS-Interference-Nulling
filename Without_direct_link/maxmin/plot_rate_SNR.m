close all
load main_SDR1.mat
power_set = -sigma2_dB_set;
l=plot(power_set, rate_set(:,1), '<-','Color',[0.00,0.45,0.74]); hold on
l.MarkerFaceColor = l.Color;
l=plot(power_set, rate_set(:,2), '<--','Color',[0.00,0.45,0.74]); hold on
l.MarkerFaceColor = l.Color;
load main_subgradient1.mat
power_set = -sigma2_dB_set;
l=plot(power_set, rate_set(:,1), 'o-','Color',[0.85,0.33,0.10]); hold on
l.MarkerFaceColor = l.Color;
l=plot(power_set, rate_set(:,2), 'o--','Color',[0.85,0.33,0.10]); hold on
l.MarkerFaceColor = l.Color;
legend('SDR+Bisection Search ($N=16$)','SDR+Bisection Search ($N=36$)',...
    'Subgradient Method ($N=16$)','Subgradient Method ($N=36$)','Location','northwest','Interpreter','latex')
xlabel('Transimit power (dBm) ','Interpreter','latex')
ylabel('Minimum rate (bps/Hz) ','Interpreter','latex')
grid on