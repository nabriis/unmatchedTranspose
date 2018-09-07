set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
figure(1)
clf
n = [64^2, 128^2, 256^2, 512^2, 1024^2, 2048^2];
k = [508, 1002, 1932, 2963, 5568, 9120];

loglog(n,k,'o-','linewidth',2.5)
hold on
loglog(n,50*n.^(1/3),'--','linewidth',2.5)
loglog(n,6*n.^(1/2),':','linewidth',2.5)
hold off
axis tight

set(gca,'TickLabelInterpreter','latex'); set(gca,'fontsize',14)
legend('Krylov--Schur','$\mathcal{O}(n^{1/3})$','$\mathcal{O}(n^{1/2})$','location','northwest')
xlabel('$n$')
ylabel('$\#$ MVM')

%%
savefig('results/ASTRAMVM.fig')
print('results/ASTRAMVM.eps','-depsc')
 