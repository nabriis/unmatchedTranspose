function y = prec(M1, M2, x, transp_flag)

%PREC  Perform action of preconditioner
% function y = prec(M1, M2, x, transp_flag)
%     y = M2  \ (M1  \ x)
% or  y = M1' \ (M2' \ x)
%
% Revision date: June 5, 2010
% (C) Michiel Hochstenbach 2016

if nargin < 4 || isempty(transp_flag)
  transp_flag = 0;
end

if ~transp_flag
  % inv(M1 M2) = inv(M2) inv(M1)
  y = minvv(M2, minvv(M1, x, 0), 0);
else
  % Transpose: inv(M1 M2)' = inv(M1') inv(M2')
  y = minvv(M1, minvv(M2, x, 1), 1);
end
