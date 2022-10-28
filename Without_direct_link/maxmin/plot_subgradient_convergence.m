close all
load subgradient_convergence.mat
plot(1:iter_max, rate_tract);
xlabel('Number of iterations', 'Interpreter','latex')
ylabel('Minimun rate (bps/Hz)', 'Interpreter','latex')