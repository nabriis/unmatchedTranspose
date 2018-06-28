function [theta, D, x, X, y, Y, KX] = extraction_ep(extraction, target, H, K, opts)

%EXTRACTION_EP  Various extraction processes for the (generalized) eigenvalue problem
% function [theta, D, x, X, y, Y, KX] = extraction_ep(extraction, target, H, K, opts)
%
% Input arguments:
%   extraction : standard | harmonic | refined | relative | rational | minres | regularized
%   target     : complex number | '-inf' | 'inf' | 'real' | 'imag' | 'abs'
%   H          : U'*A*U  (V'*A*U for two-sided)
%   K          : U'*B*U  (V'*B*U for two-sided)
% Other arguments that may be needed:
%   fix        : if fix then base extraction on target otherwise on Rayleigh quotient
%   nr         : number of vectors required
%   a, b       : weights for regularized extraction
% extraction     one-sided           | two-sided
%   standard     --                  | --
%   harmonic     AtBU, AU, BU        | AtBV, ACV, BCV, VAtAU
%   refined      AtBU, AU, BU        | AtBV, ACV, BCV
%   relative     AtBU, AU, BU        | AtBV, ACV, BCV, VAtAU
%   rational     AtBU, AU, BU, AtcBU | AtBV, ACV, BCV, VAtAU
%   minres       AU, BU              | ACV, BCV
%   regularized  UACAU               | UAU, VACV, VAACV, UU, VV
% where
%   ACV        = A'*V
%   BCV        = B'*V
%   AtBU       = AU-target*BU = QR
%   AtBV       = (A-target*B)'*V
%   UACAU      = U'*A'*A*U
%   VAtAU      = V'*(A-target*B)^2*U
%   AtcBU      = AU+target'*BU
%
% Output arguments:
%   theta : approximate eigenvalue
%   D     : approximate eigenvalues
%   x     : primitive right Ritz/Petrov vector
%   X     : primitive right Ritz/Petrov vectors
%   y     : primitive left  Ritz/Petrov vector
%   Y     : primitive left  Ritz/Petrov vectors
%
% References:
%   G.W. Stewart, Matrix Algorithms, Volume II: Eigensystems, SIAM 2001
%
%   M.E. Hochstenbach, 
%     Generalizations of harmonic and refined Rayleigh-Ritz,
%     ETNA 20, pp. 235-252, 2005.
%
%   M.E. Hochstenbach, 
%     Variations on harmonic Rayleigh-Ritz for the generalized and polynomial eigenvalue problems,
%     Preprint, Case Western Reserve University, September 2005
%
% See also EXTRACTION_QEP, EXTRACTION_PEP, JD
%
% Revision date: January 22, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 4
  K = [];
end
if nargin < 5
  opts = [];
end

twosided = getfield0(opts, 'twosided', 0);
fix      = getfield0(opts, 'fix', 1);
nr       = getfield0(opts, 'nr', []);

if isa(extraction, 'function_handle')
  % maximize function
  [X, D] = eig0(H, K);
  [~, index] = sort(1./abs(extraction(D)));
  D = D(index);
  X = X(:,index);
  x = X(:,1);
  theta = D(1);
  return
end

if ~fix || strcmp(extraction, 'standard') || strcmp(extraction, 'minres') || ...
           strcmp(extraction, 'regularized')
  if ~twosided
    [X, S] = qz0(H, K, target);
    D = diag(S);
    % [X, D] = eig0(H, K, target);
    % [X, D] = eig_spd(H, K, target); % for optimal stability if K SPD
  else
    [X, D, Y, KX] = eig0(H, K, target);
  end
end

if strcmp(extraction, 'harmonic')
  if ~twosided
    if ~fix % use Ritz value
      [Q,~] = qr_pos(opts.AU-D(1)*opts.BU, 0);
      [X,D] = qz0(Q'*opts.BU, opts.R, 'abs');
    else
      [X,D] = qz0(opts.QBU, opts.R, 'abs');
      % [X,D] = eig0(opts.QBU, opts.R, 'abs');
    end
  else
    if ~fix               % Use Ritz value
      AtBU  = opts.AU -D(1) *opts.BU;
      AtBV  = opts.ACV-D(1)'*opts.BCV;
      VAtAU = AtBV'*AtBU;
      if isempty(K)
        HtK = H-D(1)*eye(size(H));
      else
        HtK = H-D(1)*K;
      end
      [X, D, Y, KX] = eig0(HtK, VAtAU, 'abs');
    else                  % Use target
      if isempty(K)
        HtK = H-target*eye(size(H));
      else
        HtK = H-target*K;
      end
      [X, D, Y, KX] = eig0(HtK, opts.VAtAU, 'abs');
    end
  end
end

if strcmp(extraction, 'relative')
  if ~twosided
    if ~fix                % Use Ritz value
      [Q,~] = qr_pos(opts.AU-D(1)*opts.BU, 0);
      [X,D] = qz0(Q'*opts.AU, opts.R, 'abs');
    end
      [X,D] = qz0(opts.QAU, opts.R, 'abs');
      % [X,D] = eig0(opts.QAU, opts.R, 'abs');
  else
    if ~fix               % Use Ritz value
      AtBU  = opts.AU -D(1) *opts.BU;
      AtBV  = opts.ACV-D(1)'*opts.BCV;
      VAtAU = AtBV'*AtBU;
      [X,D,Y] = eig0(AtBV'*opts.AU, VAtAU, 'abs');
    end
    [X,D,Y] = eig0(opts.AtBV'*opts.AU, opts.VAtAU, 'abs');
  end
end

if strcmp(extraction, 'rational')
  if ~twosided
    if ~fix % Use Ritz value
      [Q,~] = qr_pos(opts.AU-D(1)*opts.BU, 0);
      [X,D] = qz0(Q'*opts.AU+target'*(Q'*opts.BU), opts.R, 'abs');
    end
    [X,D] = qz0(opts.QAU+target'*opts.QBU, opts.R, 'abs');
    % [X,D] = eig0(opts.QAU+target'*opts.QBU, opts.R, 'abs');
  else
    if ~fix               % Use Ritz value
      AtBU = AU -D(1) *BU;
      AtBV = ACV-D(1)'*BV;
      VAtAU = AtBV'*AtBU;
      [X,D,Y] = eig0(AtBV'*(AU+target'*BU), VAtAU, 'abs');
    end
    [X,D,Y] = eig0(opts.AtBV'*(opts.AU+target'*opts.BU), opts.VAtAU, 'abs');
  end
end

if strcmp(extraction, 'refined')
  if ~fix || ischar(target)   % Use Ritz value
    [D, X] = svd0(opts.AU-D(1)*opts.BU, 0, 0, 1);
  else
    [D, X] = svd0(opts.R, 0, 0);
  end
  if twosided                 % Only possible with orthonormal bases
    if ~fix || ischar(target) % Use Ritz value
      [D2, Y] = svd0(opts.ACV-D(1)'*opts.BCV, 0, 0, 1);
    else
      [D2, Y] = svd0(opts.AtBV, 0, 0, 1);
    end
  end
end

if strcmp(extraction, 'minres')
  if ~twosided
    [theta, x] = extraction_minres(opts.AU, opts.BU, H, [], D(1), [], 1e-1);
  else
    [theta, x, y] = extraction_minres(opts.AU, opts.BU, H, [], D(1), [], 1e-1, opts.ACV, opts.BCV);
  end
end

if strcmp(extraction, 'regularized')
  opts.a = getfield0(opts, 'a', 1);
  opts.b = getfield0(opts, 'b', 1);
  if ~twosided
    [theta, x] = extraction_ep_reg(target, H, [], opts.UACAU, D(1), X(:,1), opts.a, opts.b, [], 1e-2);
  else
  	[theta, x, y] = extraction_ep_reg(target, H, K, opts.UACAU, D(1), X(:,1), opts.a, opts.b, [], 1e-2, ...
      opts.VAACV, opts.UAU, opts.VACV, opts.UU, opts.VV, Y(:,1));
  end
end

if ~isempty(nr)
  X = X(:,1:nr);
  D = D(:,1:nr);
  if twosided
    Y = Y(:,1:nr);
  end
end

if ~strcmp(extraction, 'minres') && ~strcmp(extraction, 'regularized')
  x = X(:,1);
  if twosided
    y = Y(:,1);
  end
end

if strcmp(extraction, 'standard')
  theta = D(1);
else                       % Take extra Rayleigh quotient
  if ~twosided
    theta = rq(H,K,x);
  else
    theta = grq(H,K,x,y);
  end
end

% Ensure realness if applicable
if strcmp(target, 'infreal') || strcmp(target, '-infreal')
  [theta, D, x, X] = eig_real(theta, D, x, X);
end
