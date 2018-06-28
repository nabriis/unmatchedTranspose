function [x, flag] = check_relres(x, y, relres)

%CHECK_RELRES  If relres of iterative method >= 1 then take default vector (usually right-hand side)
% function [x, flag] = check_relres(x, y, relres)
%
% See also CHECK_FLAG
%
% Revision date: February 20, 2006
% (C) Michiel Hochstenbach 2002-2006

flag = 0;
if relres >= 1
  if nargout < 2
    fprintf('Info: relres = %g, taking residual instead\n', relres);
  else
    flag = 1;
  end
  x = y;
end
