function [c, flag] = divide(a, b, str, tol)

%DIVIDE  Division without any checks
% function [c, flag] = divide(a, b, str, tol)
%
% See also DIVIDE_NEAT
%
% Revision date: March 26, 2004
% (C) Michiel Hochstenbach 2002-2004

flag = 0;
c = a/b;
