clc; close all; clear;
%%
load('results/PedagogicalExample.mat')
%%
%prepare for plotting
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%%
%temp stuff
%title(sprintf('$ \\min \\Re(\\lambda) = %.8f, \\max \\Re(\\lambda) = %.6f $',min(real(eigenvalues_sr)),max(real(eigenvalues_lr))))
%%
%Figure 1 (BA w.r.t. xn_ast)
figure(1);
semilogy(kk,err_xn_BA_well_xn_ast,'linewidth',3)
hold on
semilogy(kk,err_xn_BA_ill_xn_ast,'--','linewidth',3)
hold off
set(gca,'TickLabelInterpreter','latex'); set(gca,'fontsize',16)
ylim([10^(-16) 10^(16)])
xticks([10^5 3*10^5, 5*10^5, 7*10^5 9*10^5])
yticks([10^(-12) 10^0 10^12])
legend('$A_{\mathrm{well}}$','$A_{\mathrm{ill}}$','location','northwest')
xlabel('Iteration number $k$')
ylabel('$ \| x^k - x^\ast\| / \|x^\ast\|$','fontsize',20)

%%
%Figure 2 (SBA w.r.t. xn_alpha_ast)
figure(2);
semilogy(kk,err_xn_SBA_well_xn_alpha_ast,'linewidth',3)
hold on
semilogy(kk,err_xn_SBA_ill_xn_alpha_ast,'--','linewidth',3)
hold off
set(gca,'TickLabelInterpreter','latex'); set(gca,'fontsize',16)
ylim([10^(-16) 10^(16)])
xticks([10^5 3*10^5, 5*10^5, 7*10^5 9*10^5])
yticks([10^(-12) 10^0 10^12])
legend('$A_{\mathrm{well}}$','$A_{\mathrm{ill}}$','location','northwest')
xlabel('Iteration number $k$')
ylabel('$ \| x^k - x_\alpha^\ast\| / \|x_\alpha^\ast\|$','fontsize',20)


%%
%Figure 3 (BA and SBA w.r.t xex for well-cond)
figure(3);
loglog(kk,err_x_BA_well_xex,'linewidth',3)
hold on
loglog(kk,err_xn_BA_well_xex,':','linewidth',3)
loglog(kk,err_x_SBA_well_xex,'--','linewidth',3)
loglog(kk,err_xn_SBA_well_xex,'-.','linewidth',3)
hold off
set(gca,'TickLabelInterpreter','latex'); set(gca,'fontsize',16)
ylim([10^(-4) 10^(3)])
xticks([10^0, 10^2, 10^4])
legend('BA $e=0$','BA $e\not=0$','SBA $e=0$','SBA $e\not=0$','location','northwest')
xlabel('Iteration number $k$')
ylabel('$ \| x^k - \bar{x}\| / \|\bar{x}\|$','fontsize',20)

%%
%Figure 4 (BA and SBA w.r.t xex for ill-cond)
figure(4);
loglog(kk,err_x_BA_ill_xex,'linewidth',3)
hold on
loglog(kk,err_xn_BA_ill_xex,':','linewidth',3)
loglog(kk,err_x_SBA_ill_xex,'--','linewidth',3)
loglog(kk,err_xn_SBA_ill_xex,'-.','linewidth',3)
hold off
set(gca,'TickLabelInterpreter','latex'); set(gca,'fontsize',16)
ylim([10^(-4) 10^(3)])
xticks([10^0, 10^2, 10^4])
legend('BA $e=0$','BA $e\not=0$','SBA $e=0$','SBA $e\not=0$','location','northwest')
xlabel('Iteration number $k$')
ylabel('$ \| x^k - \bar{x}\| / \|\bar{x}\|$','fontsize',20)

%%
figure(5);
loglog(kk,err_x_BA_well_x_ast,'b','linewidth',3)
hold on
loglog(kk,err_x_SBA_well_x_alpha_ast,'k','linewidth',3)
loglog(kk,err_x_BA_well_xex,'b--','linewidth',3)
loglog(kk,err_x_SBA_well_xex,'k--','linewidth',3)

%%
figure(6);
loglog(kk,err_xn_BA_well_xn_ast,'b','linewidth',3)
hold on
loglog(kk,err_xn_SBA_well_xn_alpha_ast,'k','linewidth',3)
loglog(kk,err_xn_BA_well_xex,'b--','linewidth',3)
loglog(kk,err_xn_SBA_well_xex,'k--','linewidth',3)
%%
% figure(1)
% savefig('results/BAe.fig')
% print('results/BAe.eps','-depsc')
% 
% figure(2)
% savefig('results/SBAe.fig')
% print('results/SBAe.eps','-depsc')
% 
% figure(3)
% savefig('results/SemiWell.fig')
% print('results/SemiWell.eps','-depsc')
% 
% figure(4)
% savefig('results/SemiIll.fig')
% print('results/SemiIll.eps','-depsc')


%%
format short e
norm(x_ast_well - xex) 
norm(x_alpha_ast_well - xex) 
norm(x_alpha_ast_well - x_ast_well) 

norm(x_ast_ill-xex) 
norm(x_alpha_ast_ill-xex) 
norm(x_alpha_ast_ill-x_ast_ill) 

norm(xn_ast_well - xex) 
norm(xn_alpha_ast_well - xex) 
norm(xn_alpha_ast_well - xn_ast_well) 

norm(xn_ast_ill-xex) 
norm(xn_alpha_ast_ill-xex) 
norm(xn_alpha_ast_ill-xn_ast_ill) 