function rho = rq(A, B, u, normalize)

%RQ  Rayleigh Quotient of A, B and u
% function rho = rq(A, B, u, normalize)
%   rho = u'Au/u'Bu
% when called with 2 arguments, it computes rho = u'Au/u'u
%
% See also GRQ, PRQ
%
% Revision date: January 22, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 4 || isempty(normalize)
  normalize = 1;  % Normalization necessary
end

if nargin == 2           % Rayleigh quotient of A and u
  if normalize
    B = B / norm(B);
  end
  rho = B'*mv(A, B);
else                     % generalized eigenvalue problem
  rho = (u'*mv(A, u)) / (u'*mv(B, u));
end
