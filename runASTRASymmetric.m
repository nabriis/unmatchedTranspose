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

%Set-up ASTRA projector, solution and data
A = paralleltomo_astra(N,linspace(180/t,180,t),p,N,'cuda');
B = A'; [m,n] = size(A);

A = A*eye(n,n); %Generate matrix
B = B*eye(m,m); %Generate matrix

%%
BA = B*A;

%%
norm(BA*BA'-BA'*BA,'fro')/norm(BA,'fro')^2

%%
%BA_sym = 0.5*((BA)'+(BA));
BA_skew = 0.5*((BA)'-(BA));
%%
norm(BA_skew,'fro')/norm(BA,'fro')

