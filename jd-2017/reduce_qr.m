function [Q, R, Q1] = reduce_qr(Q, R, X)

%REDUCE_QR  Reduce QR
% function [Q, R, Q1] = reduce_qr(Q, R, X)
%
% Revision date: February 3, 2006
% (C) Michiel Hochstenbach 2014

% instead of [Q, R] = qr_pos(Q*R*X)
[Q1, R] = qr_pos(R*X, 0);
Q = Q*Q1;
