%clc; 
clear; close all force;
% Estimate the left-most eigenvalue (real part)
% of BA for GPU ASTRA projectors.
rng(0)

%Add paths
addpath(genpath('astra-1.8'))
addpath(genpath('spot-master'))
addpath(genpath('Functions'))
addpath(genpath('ks'))
addpath(genpath('jd-2017'))

format short e

%Paramters for tests
matrix = 0; %Use matrix representation?, 0=no, 1=yes.
n_trails = 5; %Number of trails to average over.
method = 'ks'; %eigs, ks, jd, symm.
proj_type = 'cuda'; %GPU 'cuda', change to 'line' for CPU version

%Other
maxit = 1500;
tol = 10^(-2);
maxit_symm = 5; %Max iterations for symmetric method.
tol_symm = 10^(-6);

% -------- Parameters for CT Problem ---------
%Parameters defining matrix size (n=N*N, m=t*p)
%N: NxN image domain,
%t: Number of measurement locations
%p: Number of detector pixels per measurement
%N=64; t = 42; p = 39; %-0.6    %Symm method is not better for any case.
%N=64; t = 70; p = 42; %-0.11   %Try symm method here w. 10 max_it.
N=128; t = 90; p = 80; %-0.928 %Symm method is not better for any case.
%N=128; t = 140; p = 84; %-0.112 %Try symm w. 15, or 25 max_it
%N=256; t = 200; p = 160; %-1.640 %Symm method is not better for any case.
%N=512; t = 360; p = 320;
%N=1024; t = 720; p = 640;
%N=2048; t = 1440; p = 1280;

%CT specific parameters (should be no need to change)
theta = linspace(180/t,180,t); %Measurement locations in degrees
domain = N; dl = N;

%Generate matrix or matrix representation (ASTRA)
disp('Generating System..')
if matrix == 0
    A = paralleltomo_astra(N,theta,p,dl,proj_type);
    B = A'; [m,n] = size(A);
elseif matrix==1
    A = paralleltomo_astra(N,theta,p,dl,proj_type);
    B = A'; [m,n] = size(A);
    A = A*eye(n,n); %Generate matrix
    B = B*eye(m,m); %Generate matrix
end

%Method parameters
opts.maxit = maxit; opts.tol = tol;
opts.mindim = 30; opts.maxdim = 60;
opts.nr = 1; opts.target = '-inf';
opts.n = size(A,1); opts.m = size(A,1);
opts.absrel = 'abs';

%jd sepcific
opts_jd = opts;
%opts_jd.solver = 'bicgstab';
%opts_jd.maxit_inner = 5;
opts_jd.krylov = 1; opts_jd.thick=1;

%Symm specific
opts_symm = opts;
opts_symm.maxit = maxit_symm;
opts_symm.tol = tol_symm;

eigenvalue_sr = zeros(n_trails,1);


disp('Running trails..')

for i = 1:n_trails
    
    %Timing starts
    tic;
    
    switch method
        case 'eigs'
            [eigenvalues_eigs] = eigs(@(x) A*(B*x),m,1,'sr',opts);
            eigenvalue_sr(i) = min(real(eigenvalues_eigs));
            n_BA(i) = Inf;
        case 'ks'
            [eigenvalues_ks, ~, ~,~,nr_mv] = krylov_schur(@(x) A*(B*x), opts);
            n_BA(i) = 2*nr_mv(end);
            eigenvalue_sr(i) = min(real(eigenvalues_ks(1,1)));
        case 'jd'
            [eigenvalues_jd,~,~,nr_mv] = jd(@(x) A*(B*x), [], opts_jd);
            n_BA(i) =2*nr_mv(end);
            eigenvalue_sr(i) = min(real(eigenvalues_jd));
        case 'symm'
            [H,~,~,~,nr_mv] = krylov_schur(@(x) A*(B*x), opts_symm);
            Hs = H(1:end-1,:);
            n_BA(i) = 2*nr_mv(end);
            eigenvalues_symm = eig(0.5*(Hs+Hs'));
            eigenvalue_sr(i) = H(1,1);
            eigenvalue_sr2(i) = min(real(eigenvalues_symm));
            
    end
    
    t_n(i) = toc;
    %Timing ends
    
end
%%
disp(' ')
fprintf('Timing for %s.. mean: %3.4f seconds, std: %3.4f seconds\n',method,mean(t_n),std(t_n));
try
fprintf('Iterations for %s.. mean: %3.4f, std: %3.4f\n',method,mean(n_BA),std(n_BA));
catch
end
if strcmp(method,'symm')
    fprintf('Left-most eigenvalue using real(H(1,1))                 mean: %6.9f, std: %6.9f\n',mean(eigenvalue_sr),std(eigenvalue_sr));
    fprintf("Left-most eigenvalue using min(real(eig(0.5*(Hs+Hs')))) mean: %6.9f, std: %6.9f\n",mean(eigenvalue_sr2),std(eigenvalue_sr2));
else
    fprintf('Left-most eigenvalue mean: %6.9f, std: %6.9f\n',mean(eigenvalue_sr),std(eigenvalue_sr));
end
disp(' ')

