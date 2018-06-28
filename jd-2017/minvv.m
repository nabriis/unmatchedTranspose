function y = minvv(A, x, transp_flag)

%MINVV  Matrix inverse times vector: y = A\x
% function y = minvv(A, x, transp_flag)
%
% Revision date: January 22, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 3 || isempty(transp_flag)
  transp_flag = 0;
end

if isnumeric(A)
  if ~isempty(A)
    if ~transp_flag
      y = A\x;
    else
      y = A'\x;
    end
  else
    % A = I
    y = x;
  end
else
  % Assume inv(A) performed by function A
  y = mv(A, x, transp_flag);
end
