function x = mv(A, x, transp_flag)

%MV  Compute matrix times vector
% function x = mv(A, x, transp_flag)
%
% This function is used in iterative routines to ensure that functions
%   carrying out a matrix-vector product can be used
%
%   A   matrix or function, empty matrix is seen as the identity
%   x   vector or matrix with matching dimension
%
% Revision date: January 15, 2014
% (C) Michiel Hochstenbach 2016

if nargin < 3 || isempty(transp_flag)
  transp_flag = 0;
end
if isnumeric(A) || ~isa(A,'function_handle')
  if ~isempty(A)  % Empty matrix A stands for identity
    if ~transp_flag
      x = A*x;
    else 
      x = A'*x;   % (x'*A)' used to be faster; not any longer
    end
  end
else              % Function
  if nargin < 3
    x = A(x);
  else
    if nargin(A) == 1
      if transp_flag
        warning('Transpose flag set for function with 1 argument');
      else
        x = A(x);
      end
    else
      x = A(x, transp_flag);
    end
  end
end
