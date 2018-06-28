function x = projss(x, U, V, VU)

%PROJSS  Orthogonal or oblique projection onto a subspace
% function x = projss(x, U, V, VU)
%
% y = (I-U (V'U)\U') x, projection onto V_perp along U
% Assumes V = U if not given: orthogonal projection
% Assumes VU (= V'U) = I if not given: orthonormal basis or biorthonormal bases
%
% NB: This is Gram--Schmidt, projecting once instead of twice
% This procedure is cheap but may be numerically unstable
% It may be sufficient, for instance in generating a Krylov space K((I-UU')A, b)
% However, use RGS instead if stability is necessary
%
% See also RGS
%
% Revision date: August 3, 2016
% (C) Michiel Hochstenbach 2016

if ~isempty(U) && ~isempty(x)
  if nargin == 2 || nargin == 3 && isempty(V) || nargin == 4 && isempty(V) && isempty(VU)
    x = x - U*(U'*x);        % I-UU'
  elseif nargin == 3 && size(V,1) < size(U,1)
    x = x - U*(V \ (U'*x));  % I-U(U'U)\U', here V is used for U'U
  elseif nargin == 3
    x = x - U*(V'*x);        % I-UV'
  else
    x = x - U*(VU \ (V'*x)); % I-U(V'U)\V'
  end
end
