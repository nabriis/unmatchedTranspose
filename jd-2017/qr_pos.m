function [Q,R] = qr_pos(A, full)

%QR_POS  QR with positive diagonal of R
% function [Q,R] = qr_pos(A, full)
%
% - Matlab's qr may return a [Q,R] combination with negative
%   diagonal entries (in particular R(1,1) < 0).
%   If this is the case, flip the columns.
% - Moreover, if called with one outout argument, return
%   Q instead of R (as Matlab does).
%
% See also QR
%
% Revision date: April 20, 2005
% (C) Michiel Hochstenbach 2002-2005

if nargin == 1
  [Q,R] = qr(A);
else
  [Q,R] = qr(A, 0);
end

for j = 1:min(size(R))
  if R(j,j) < 0
    Q(:,j) = -Q(:,j);
    R(j,:) = -R(j,:);
  end
end
