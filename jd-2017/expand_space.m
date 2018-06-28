function [U, H] = expand_space(U, S)

%EXPAND_SPACE  Expand an orthonormal basis with extra vector(s)
% function [U, H] = expand_space(U, S)
% In:  U'*U = I, U: n times k, S: n times l
% Out: U'*U = I
%
% See also RGS
%
% Revision date: June 10, 2007
% (C) Michiel Hochstenbach 2014

k = size(U,2);
for j = 1:size(S,2)
  if nargout > 1
    [U(:,k+j), H(1:k+j,j)] = rgs(S(:,j), U);
  else
    U(:,k+j) = rgs(S(:,j), U);
  end
end
