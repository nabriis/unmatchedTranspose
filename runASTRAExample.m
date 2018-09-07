clc; clear; close all; rng(0)

addpath(genpath('astra-1.8'))
addpath(genpath('spot-master'))
addpath(genpath('Functions'))

%Parameters defining matrix (n=N*N, m=t, pvec)
%N=64; t = 42; p = 39;
%N=64; t = 70; p = 42;
N=128; t = 90; p = 80;
%N=256; t = 200; p = 160;
%N=512; t = 360; p = 320;
%N=1024; t = 720; p = 640; % -5.67589
%N=2048; t = 1440; p = 1280; % -??

%noise level
rnl_e = 0.05;

%Set-up ASTRA projector, solution and data
A = paralleltomo_astra(N,linspace(180/t,180,t),p,N,'cuda');
B = A'; [m,n] = size(A);

%data
xex=phantom(N); xex = xex(:);
b = A*xex;
e=randn(length(b),1);
bn = b+rnl_e*norm(b)*e/norm(e);

disp('Calculating eigenvalues using eigs')
opts.maxit = 500; opts.tol = 10^(-2);
eig_lr = eigs(A*B,1,'lr',opts);
eig_sr = eigs(A*B,1,'sr',opts);

%%
%Shift eigenvalues by factor 2 min eig
alpha2 = 2*abs(min(real(eig_sr))); 

disp('Calculating limit solutions')
%Calculate solutions that we should convergece to (need to generate matrix)
A_m = A*speye(size(A,2));
B_m = B*speye(size(B,2));

x_ast = B_m*((A_m*B_m)\b); 
xn_ast = B_m*((A_m*B_m)\bn);
x_alpha_ast = B_m*((A_m*B_m+alpha2*eye(m,m))\b);
xn_alpha_ast = B_m*((A_m*B_m+alpha2*eye(m,m))\bn);

clear('A_m'); clear('B_m')

%%

%Calculate step-size (w) for BA and Reg-BA
factor = 1.9; %Step size factor, must be in (0,2)

%Step size estimate for BA iteration
wBA = factor/(normest(B*A)); %w = factor/(normest(A)^2)

%Estimate step size (!!slow!!)
disp('Estimating step size for SBA')
wSBA = est_step_size_mf(A,B,alpha2,factor);
wSBA_f = est_step_size_fast(A,B,alpha2,factor); 

%%
%parameters for iterations
max_iter= 10^6;
it=1; dk = 100; kk = unique(round(logspace(0,log10(max_iter),dk)));

%Preallocate a few things
nx = norm(xex);
nx_ast = norm(x_ast);
nxn_ast = norm(xn_ast);
nx_alpha_ast = norm(x_alpha_ast);
nxn_alpha_ast = norm(xn_alpha_ast);

%All the solutions
xBA = zeros(n,1); xnBA = zeros(n,1);
xSBA = zeros(n,1); xnSBA = zeros(n,1);
xSBA_f = zeros(n,1); xnSBA_f = zeros(n,1);
disp('Running iterations')
%Run iterations
for k=1:max_iter
    %BA
    g = B*(b-A*xBA); 
    xBA = xBA + wBA*g;
    
    g = B*(bn-A*xnBA); 
    xnBA = xnBA + wBA*g;
    
    %SBA
    g = B*(b-A*xSBA);
    xSBA = (1-alpha2*wSBA)*xSBA + wSBA*g;
    
    g = B*(bn-A*xnSBA);
    xnSBA = (1-alpha2*wSBA)*xnSBA + wSBA*g;

    %SBA fast
    g = B*(b-A*xSBA_f);
    xSBA_f = (1-alpha2*wSBA_f)*xSBA_f + wSBA_f*g;
    
    g = B*(bn-A*xnSBA_f);
    xnSBA_f = (1-alpha2*wSBA_f)*xnSBA_f + wSBA_f*g;
    
    if any( k==kk )
        k/max_iter
        
        err_x_BA_xex(it) = norm(xBA-xex)/nx;
        err_xn_BA_xex(it) = norm(xnBA-xex)/nx;
        err_x_BA_x_ast(it) = norm(xBA-x_ast)/nx_ast;
        err_xn_BA_xn_ast(it) = norm(xnBA-xn_ast)/nxn_ast;

        err_x_SBA_xex(it) =  norm(xSBA-xex)/nx;
        err_xn_SBA_xex(it) =  norm(xnSBA-xex)/nx;
        err_x_SBA_x_alpha_ast(it) =  norm(xSBA-x_alpha_ast)/nx_alpha_ast;
        err_xn_SBA_xn_alpha_ast(it) =  norm(xnSBA-xn_alpha_ast)/nxn_alpha_ast;      
        
        err_x_SBA_f_xex(it) =  norm(xSBA_f-xex)/nx;
        err_xn_SBA_f_xex(it) =  norm(xnSBA_f-xex)/nx;
        err_x_SBA_f_x_alpha_ast(it) =  norm(xSBA_f-x_alpha_ast)/nx_alpha_ast;
        err_xn_SBA_f_xn_alpha_ast(it) =  norm(xnSBA_f-xn_alpha_ast)/nxn_alpha_ast;      
        
        it = it+1;
    end
end

%%
%save('results/ASTRAExample2.mat')