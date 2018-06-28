clc; close all; clear;
%%
load('results/ASTRAExample.mat')
%%
%prepare for plotting
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%%
%temp stuff
%title(sprintf('$ \\min \\Re(\\lambda) = %.8f, \\max \\Re(\\lambda) = %.6f $',min(real(eigenvalues_sr)),max(real(eigenvalues_lr))))
%%
%Figure 1 (BA w.r.t. x_ast)
figure(1);
loglog(kk,err_x_BA_x_ast,'--','linewidth',3)
hold on
loglog(kk,err_x_SBA_x_alpha_ast,'linewidth',3)
hold off
set(gca,'TickLabelInterpreter','latex'); set(gca,'fontsize',16)
ylim([10^(-5) 10^(2)])
xlim([0 10^5])
xticks([10^0, 10^1, 10^2, 10^3, 10^4, 10^5])
h=legend('BA: \hspace{0.2cm}$\| x^k - \bar{x}^*\| / \| \bar{x}^*\|$','SBA: $\| x^k - \bar{x}_\alpha^*\| / \| \bar{x}_\alpha^*\|$','location','northwest');
set(h,'FontSize',18)
xlabel('Iteration number $k$')
%ylabel('$ \| x^k - \bar{x}\| / \|\bar{x}\|$','fontsize',20)

%%
%Figure 5 (BA and SBA w.r.t xex for well-cond)
figure(2);
loglog(kk,err_xn_BA_xex,'linewidth',3)
hold on
loglog(kk,err_x_BA_xex,':','linewidth',3)
loglog(kk,err_x_SBA_xex,'--','linewidth',3)
loglog(kk,err_xn_SBA_xex,'-.','linewidth',3)
hold off
set(gca,'TickLabelInterpreter','latex'); set(gca,'fontsize',16)
ylim([0.25 10^(0)])
xlim([0 10^5])
xticks([10^0, 10^1, 10^2, 10^3, 10^4, 10^5])
legend('BA $e=0$','BA $e\not=0$','SBA $e=0$','SBA $e\not=0$','location','north')
xlabel('Iteration number $k$')
ylabel('$ \| x^k - \bar{x}\| / \|\bar{x}\|$','fontsize',20)
%%

% figure(1)
% savefig('results/ASTRAnoisefree.fig')
% print('results/ASTRAnoisefree.eps','-depsc')
% 
% figure(2)
% savefig('results/ASTRAsemi.fig')
% print('results/ASTRAsemi.eps','-depsc')