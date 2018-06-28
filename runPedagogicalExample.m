clc; clear; close all; rng(0)

addpath(genpath('regu'))
addpath(genpath('Functions'))
addpath(genpath('results'))

n = 64; m =64; %Matrix dimensions
rnl_E = 0.05; %Relative perturbation on A'
rnl_e = 0.05; %Relative perturbatio on b

%Generate A (max eigenvalue of A'*A = 1)
s_well = logspace(0,-2,min(m,n));
s_ill = logspace(0,-4,min(m,n));

[A_well,~,~] = regutm(m,n,s_well);
[A_ill,~,~] = regutm(m,n,s_ill);

%Generate B
E = randn(n,m);
B_well = A_well'+rnl_E*norm(full(A_well'))*E/norm(E);
B_ill  = A_ill'+rnl_E*norm(full(A_well'))*E/norm(E);

%Generate xex and b
[~,~,xex] = shaw(n);
b_well = A_well*xex; b_ill = A_ill*xex; 
e=randn(length(b_well),1);
bn_well = b_well+rnl_e*norm(b_well)*e/norm(e);
bn_ill  = b_ill+rnl_e*norm(b_ill)*e/norm(e);

disp('Calculating eigenvalues using eigs')
opts.maxit = 500; opts.tol = 10^(-2);
A_well_eig_sr = eigs(A_well*B_well,1,'sr',opts);
A_well_eig_lr = eigs(A_well*B_well,1,'lr',opts);
A_ill_eig_sr  = eigs(A_ill*B_ill,1,'sr',opts);
A_ill_eig_lr  = eigs(A_ill*B_ill,1,'lr',opts);


%%
%Shift eigenvalues by factor 2 min eig
alpha2_well = 200*abs(min(real(A_well_eig_sr))); 
alpha2_ill  = 200*abs(min(real(A_ill_eig_sr))); 

%Calculate solutions that we should convergece to
x_ast_well = B_well*(A_well*B_well)^(-1)*b_well;
xn_ast_well = B_well*(A_well*B_well)^(-1)*bn_well;
x_alpha_ast_well = B_well*(A_well*B_well+alpha2_well*eye(m,m))^(-1)*b_well;
xn_alpha_ast_well = B_well*(A_well*B_well+alpha2_well*eye(m,m))^(-1)*bn_well;

x_ast_ill = B_ill*(A_ill*B_ill)^(-1)*b_ill;
xn_ast_ill = B_ill*(A_ill*B_ill)^(-1)*bn_ill;
x_alpha_ast_ill = B_ill*(A_ill*B_ill+alpha2_ill*eye(m,m))^(-1)*b_ill;
xn_alpha_ast_ill = B_ill*(A_ill*B_ill+alpha2_ill*eye(m,m))^(-1)*bn_ill;
%%

%Calculate step-size (w) for BA and Reg-BA
factor = 1.9; %Step size factor, must be in (0,2)

%Step size estimate for BA iteration
wBA_well = factor/(normest(B_well*A_well)); %w = factor/(normest(A)^2)
wBA_ill  = factor/(normest(B_ill*A_ill)); %w = factor/(normest(A)^2)

%Estimate step size (!!slow!!)
wSBA_well = est_step_size(A_well,B_well,alpha2_well,factor); 
wSBA_ill  = est_step_size(A_ill,B_ill,alpha2_ill,factor);
%%
%parameters for iterations
max_iter= 10^6;
it=1; dk = 100; kk = unique(round(logspace(0,log10(max_iter),dk)));

%Preallocate a few things
nx = norm(xex);
nx_ast_well = norm(x_ast_well);
nx_alpha_ast_well = norm(x_alpha_ast_well);
nx_ast_ill = norm(x_ast_ill);
nx_alpha_ast_ill = norm(x_alpha_ast_ill);
nxn_ast_well = norm(xn_ast_well);
nxn_alpha_ast_well = norm(xn_alpha_ast_well);
nxn_ast_ill = norm(xn_ast_ill);
nxn_alpha_ast_ill = norm(xn_alpha_ast_ill);

%All the solutions
xBA_well = zeros(n,1);
xSBA_well = zeros(n,1);
xBA_ill = zeros(n,1);
xSBA_ill = zeros(n,1);

xnBA_well = zeros(n,1);
xnSBA_well = zeros(n,1);
xnBA_ill = zeros(n,1);
xnSBA_ill = zeros(n,1);

disp('Running iterations')
%Run iterations
for k=1:max_iter
    %Well conditioned
    g = B_well*(b_well-A_well*xBA_well); 
    xBA_well = xBA_well + wBA_well*g;
    
    g = B_well*(bn_well-A_well*xnBA_well); 
    xnBA_well = xnBA_well + wBA_well*g;
    
    g = B_well*(b_well-A_well*xSBA_well);
    xSBA_well = (1-alpha2_well*wSBA_well)*xSBA_well + wSBA_well*g;
    
    g = B_well*(bn_well-A_well*xnSBA_well);
    xnSBA_well = (1-alpha2_well*wSBA_well)*xnSBA_well + wSBA_well*g;
    
    
    %Ill conditioned
    g = B_ill*(b_ill-A_ill*xBA_ill); 
    xBA_ill = xBA_ill + wBA_ill*g;
    
    g = B_ill*(bn_ill-A_ill*xnBA_ill); 
    xnBA_ill = xnBA_ill + wBA_ill*g;
    
    g = B_ill*(b_ill-A_ill*xSBA_ill);
    xSBA_ill = (1-alpha2_ill*wSBA_ill)*xSBA_ill + wSBA_ill*g;
    
    g = B_ill*(bn_ill-A_ill*xnSBA_ill);
    xnSBA_ill = (1-alpha2_ill*wSBA_ill)*xnSBA_ill + wSBA_ill*g;    
    
    
    if any( k==kk )
       
        err_x_BA_well_xex(it) = norm(xBA_well-xex)/nx;
        err_xn_BA_well_xex(it) = norm(xnBA_well-xex)/nx;
        err_x_BA_ill_xex(it) = norm(xBA_ill-xex)/nx;
        err_xn_BA_ill_xex(it) = norm(xnBA_ill-xex)/nx;

        err_x_SBA_well_xex(it) =  norm(xSBA_well-xex)/nx;
        err_xn_SBA_well_xex(it) =  norm(xnSBA_well-xex)/nx;
        err_x_SBA_ill_xex(it) =  norm(xSBA_ill-xex)/nx;         
        err_xn_SBA_ill_xex(it) =  norm(xnSBA_ill-xex)/nx;
        
        err_x_BA_well_x_ast(it) = norm(xBA_well-x_ast_well)/nx_ast_well;
        err_xn_BA_well_xn_ast(it) = norm(xnBA_well-xn_ast_well)/nxn_ast_well;
        err_x_BA_ill_x_ast(it) = norm(xBA_ill-x_ast_ill)/nx_ast_ill;
        err_xn_BA_ill_xn_ast(it) = norm(xnBA_ill-xn_ast_ill)/nxn_ast_ill;
        
        err_x_SBA_well_x_alpha_ast(it) =  norm(xSBA_well-x_alpha_ast_well)/nx_alpha_ast_well;
        err_xn_SBA_well_xn_alpha_ast(it) =  norm(xnSBA_well-xn_alpha_ast_well)/nxn_alpha_ast_well;
        err_x_SBA_ill_x_alpha_ast(it) =  norm(xSBA_ill-x_alpha_ast_ill)/nx_alpha_ast_ill;         
        err_xn_SBA_ill_xn_alpha_ast(it) =  norm(xnSBA_ill-xn_alpha_ast_ill)/nxn_alpha_ast_ill;        
        
        it = it+1;
    end
end

%%
%save('results/PedagogicalExample.mat')