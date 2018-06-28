function [lambda, Ufound, normr, nr_mv, nr_prec] = jd(A, B, opts)

%JD  Jacobi-Davidson
% function [lambda, Ufound, normr, nr_mv, nr_prec] = jd(A, B, opts)
%
% Input arguments:                                                   Default
%   A            square matrix, real or complex (may also be function)
%   B            idem                                                []
% Opts can have the following fields:
%   n            = size(A) = size(B) when A or B is a function       size(A)
%   nr           number of desired eigenpairs                        1   
%   u1           initial space (may be more-dim.)                    rand(n,1)
%   tol          tolerance of the outer iteration                    1e-6
%   absrel       absolute or relative tolerance outer iteration      'rel'
%                  relative tolerance: ||AU-UL|| < tol * ||A||_1
%   mindim       minimum dimension of subspaces                      30
%   maxdim       maximum dimension of subspaces                      60
%   maxit        maximum number of outer iterations                  5000 
%   maxit_inner  maximum number of inner iterations                  10  
%   inner_tol    tolerance of the inner iteration                    0   
%   target       complex number|'-inf'|'inf'|'real'|'imag'|'abs'     'abs'
%   extraction   standard | harmonic | refined | relative | rational | minres | regularized
%                               (if target is a string: 'standard')  'refined'       
%   fix          nonnegative real number | 'inf'                     'inf'
%                fix target until residual norm < fix
%                the correction equation is of the form
%                  (I-Buu')(A-t1*B)(I-uu'B)s = -(A-t2*I)u
%                if residual norm < fix
%                  then t1 = Rayleigh quotient
%                  else t1 = target
%   M1, M2       M=M1*M2: preconditioner in factored form            [], []
%                  M1=M2=[] means no preconditioning
%   solver       gmres | fom                                         'gmres'
%   expansion    jd|olsen|davidson|rqi|augm|riccati|resit            'jd'
%   Mtype        left | right | implicit preconditioning             'left'
%   thick        thick restart                                       1
%   krylov       start method with Krylov space of dimension mindim  1
%   verbosity    output desired (0..2)                               0
%
% Output arguments:
%   lambda       found (approximate) eigenvalues
%   Ufound       found (approximate) eigenvectors
%   normr        norm of the residual (convergence history)
%   nr_mv        number of matrix-vector products
%   nr_prec      number of actions with a preconditioner
%
% References:
%   G.L.G. Sleijpen and H.A. van der Vorst, 
%     A Jacobi-Davidson iteration method for linear eigenvalue problems, 
%     SIAM J. Matrix Anal. Appl. 17:401-425, 1996
%
%   D.R. Fokkema, G.L.G. Sleijpen, and H.A. van der Vorst, 
%     Jacobi-Davidson style QR and QZ algorithms for the reduction of matrix pencils, 
%     SIAM J. Sci. Comp. 20:94-125, 1998
%
% Revision date: December 2, 2014
% (C) Michiel Hochstenbach 2014

% global ev;

extrapolation = 0;

% Set default parameters
narginchk(1,3);
if nargin < 2
  B = [];
end
if nargin < 3
  opts = [];
end
if isnumeric(A) || ~isa(A,'function_handle')
  [m,n] = size(A);
  if m ~= n, error('Matrix must be square'); end
  opts.m = n;
else
    n = opts.n;
    m = n;
end
[nr_ev, u1, tol, mindim, maxdim, maxit, maxit_inner, inner_tol, solver, expansion, ...
  target, extraction, fix, M1, M2, Mtype, thick, kryl, verbosity, absrel] = ...
  read_options(opts, 'nr', 'u1', 'tol', 'mindim', 'maxdim', 'maxit', ...
  'maxit_inner', 'inner_tol', 'solver', 'expansion', 'target', 'extraction', ...
  'fix', 'M1', 'M2', 'Mtype', 'thick', 'krylov', 'verbosity', 'absrel');
if strcmp(absrel, 'rel')
  tol = tol * norm(A,1);
end
swap = 0;
if ~isempty(B) && strcmp0(extraction, 'harmonic') && strcmp(target, 'abs')
  if verbosity
    disp(' ')
    disp('*********************************************************************************************')
    disp('*  Note: taking the inverse problem Bx = lambda Ax with harmonic extraction and target = 0  *');
    disp('*********************************************************************************************')
  end
  swap = 1;
  target = 0;
  [A, B] = change_order(A, B);
end

% Display info
if verbosity
  disp_info('Jacobi-Davidson, programmed by Michiel Hochstenbach', ...
    m, n, target, nr_ev, extraction, expansion, maxit, tol, fix, mindim, maxdim, ...
    maxit_inner, inner_tol, solver, thick, M1, Mtype, kryl);
else
  warning('off', 'MATLAB:rankDeficientMatrix');
end

% Init
s       = u1;
U       = [];
AU      = [];
BU      = [];
H       = []; % U'*A*U
K       = []; % U'*B*U
AtBU    = []; % (A-target*B)U
UACAU   = []; % AU'*AU
Q       = []; % QR factors of (A-tB)U
R       = [];
QAU     = []; % Q'AU
QBU     = []; % Q'BU
VAtAU   = []; % V'(A-tB)^2 U
lambda  = []; % Found eigenvalues
Ufound  = [];
MUfound = []; % M^{-1} Ufound
nr_it   = 0;
k       = 0;  % Dimension subspace
relres  = 0;
nthick  = 0;
AtBU_needed  = (strcmp0(extraction, 'harmonic') || strcmp0(extraction, 'refined') || ...
  strcmp0(extraction, 'relative') || strcmp0(extraction, 'rational') || strcmp0(extraction, 'preconditioned')) ...
  && isnumeric(target);
QAU_needed   = strcmp0(extraction, 'relative') || strcmp0(extraction, 'rational');
QBU_needed   = strcmp0(extraction, 'harmonic') || strcmp0(extraction, 'rational');
UACAU_needed = strcmp0(extraction, 'regularized');
small_res = ischar(target);
if verbosity > 1
  disp('iter  dim(U)       theta               normr    relres');
  disp('------------------------------------------------------');
end

% Preallocation
% if nargout > 2
%   normr = NaN(maxit, 1);
% end
% if nargout > 3
%   nr_mv = NaN(maxit, 1);
% end
% if nargout > 4
%   nr_prec = NaN(maxit, 1);
% end

while nr_it < maxit
  nr_it = nr_it + 1;
  if kryl % Start with Krylov spaces of dimension maxdim
    dk = maxdim;
    if isempty(B)
      [U, H] = krylov(A, u1, maxdim, 1);

      AU = U*H;
      U  = U(:,1:maxdim);
      H  = H(1:maxdim,:);
      if AtBU_needed
        AtBU = AU-target*U;
      end
      if UACAU_needed
        UACAU = AU'*AU;
      end
    else
      if isnumeric(target)
        if AtBU_needed
          [U, AU, BU, AtBU] = krylov_gep(A, B, target, u1, maxdim);
        else
          [U, AU, BU] = krylov_gep(A, B, target, u1, maxdim);
        end
        H = U'*AU;
      else
        [U, H, AU] = krylov(A, u1, maxdim, 0);
        BU = mv(B, U);
      end
      K = U'*BU;
    end
    if AtBU_needed
      [Q,R] = qr_pos(AtBU, 0);
      if QAU_needed
        QAU = Q'*AU;
      end
      if QBU_needed
        if isempty(B)
          QBU = Q'*U;
        else
          QBU = Q'*BU;
        end
      end
    end
    kryl = 0;
    if nargout > 3
      nr_mv(nr_it)   = (1 + ~isempty(B)) * maxdim;
    end
    if nargout > 4
      nr_prec(nr_it) = 0;
    end
  else % Expand spaces
    dk = size(s,2);
    U = expand_space(U, s);
    AU(:,k+1:k+dk) = mv_defl_skew(A, U(:,k+1:k+dk), Ufound, Ufound, Ufound, Ufound);
    if ~isempty(B)
      BU(:,k+1:k+dk) = mv_defl_skew(B, U(:,k+1:k+dk), Ufound, Ufound, Ufound, Ufound);
    end
    if nargout > 3
      if nr_it > 1
        nr_mv(nr_it) = nr_mv(nr_it-1) + 1 + ~isempty(B);
      else
        nr_mv(nr_it) = dk * (1 + ~isempty(B));
      end
    end
    if nargout > 4
      if nr_it > 1
        nr_prec(nr_it) = nr_prec(nr_it-1);
      else
        nr_prec(nr_it) = 0;
      end
    end
    H = add_row_column(U, AU, H);            % H = U'*A*U
    if ~isempty(B)
      K = add_row_column(U, BU, K);          % K = U'*B*U
    end
    if AtBU_needed
      if isempty(B)                          % AtBU = (A-target*B)*U
        AtBu = AU(:,k+1:k+dk)-target*U(:,k+1:k+dk);
      else
        AtBu = AU(:,k+1:k+dk)-target*BU(:,k+1:k+dk);
      end
      AtBU(:,k+1:k+dk) = AtBu;               % AtBU = (A-target*B)*U
      [Q(:,1:k+dk), R(1:k+dk,k+1:k+dk)] = expand_space(Q, AtBu);
      if QAU_needed
        QAU = add_row_column(Q, AU, QAU);    % Q'*A*U
      end
      if QBU_needed
        if isempty(B)
          QBU = add_row_column(Q, U, QBU);   % Q'*B*U
        else
          QBU = add_row_column(Q, BU, QBU);
        end
      end
    end
    if UACAU_needed
      UACAU = add_row_column(AU, AU, UACAU); % U'*A'*A*U
    end
  end
  k = k + dk;
  
  % Extraction
  if nr_it > 1
    if k > size(X,1)
      nthick = min([thick k-1]);
    else
      nthick = 0;
    end
    xold = X(:,1:nthick);
  end
  extr_opts.fix = ~small_res;
  extr_opts = write_options(extr_opts, Q, R, QAU, QBU, UACAU);
  if ~extr_opts.fix || strcmp0(extraction, 'minres')
    extr_opts.AU = AU;
    if isempty(B)
      extr_opts.BU = U;
    else
      extr_opts.BU = BU;
    end
  end
%   if strcmp0(extraction, 'standard') && k == maxdim
%     extraction2 = 'harmonic'
%   end
  if strcmp0(extraction, 'regularized')
    extr_opts.a = 1;
%     if nr_it == 1
      extr_opts.b = 1;
%     else
%       extr_opts.b = normr(nr_it-1) / normA1;
%     end
  end

  [theta, D, x, X] = extraction_ep(extraction, target, H, K, extr_opts);
    %keyboard
  %lambda = [lambda; theta];
  % Compute residual
  u  = U*x;
%   [subspace(U,ev)  subspace(u,ev)]
  Au = AU*x;
  if isempty(B)
    Bu = u;
    uBu = [];
  else
    Bu = BU*x;
    uBu = u'*Bu;
  end
  
  % Minimum polynomial extrapolation
%   if nr_it > 1 && ~extrapolation && normres < 10 * tol
%     extrapolation = 1;
%     uu = [];
%     dha
%   end
%   if extrapolation
%     if size(uu,2) == 6
%       dhi
%       D = diff([uu u], 1, 2);
%       u = u - uu * (D \ (u - uu(:,end))); 
%       u = u / norm(u);
%       uu = [];
%       residual = Au-theta*Bu;
%       normres1 = norm(residual)
%       Au = mv(A,u);
% %       Bu = mv(B,u);
%       theta = u'*Au;
%       residual = Au-theta*Bu;
%       normres2 = norm(residual)
%       if normres2 < normres1
%         uu = u;
%       end
%     else
%       uu = [uu u];
%     end
%   end

  residual = Au-theta*Bu;
  normres = norm(residual);
  if nargout > 2
    normr(nr_it) = normres;
  end
  %keyboard
  if (k == maxdim) || (normres < tol)
    X = qr_pos(X,0);
    %keyboard
    %nr_it+mindim
    if normres < tol || (nr_it+mindim == (maxit+1))
        if nr_it+mindim == maxit+1
            warning('Warning, reached maxit')
        end
      % found
      if swap
        theta = 1/theta;
      end
      %keyboard
      lambda = [lambda; theta];
      Ufound = [Ufound u];
      if verbosity
        fprintf('**** Nr %d found after %d iterations with norm residual = %g\n', ...
          length(lambda), nr_it, normres);
        fprintf('lambda = %9g %+9gi\n', real(theta), imag(theta));
      end
      %keyboard
      if length(lambda) == nr_ev %|| abs(imag(theta)) < 1e-7
        %keyboard
        return
      end
      found = 1;
      MUfound = [MUfound prec(M1, M2, u)];
      if nargout > 4
        nr_prec(nr_it) = nr_prec(nr_it) + ~isempty(M1);
      end
      % new best vector
      if strcmp0(extraction, 'standard')
        x = X(:,2);
        theta = D(2);
      else
        % x1'*x2 = 0
        % x = proj(X(:,2), x);
        x = X(:,2);
        % take extra Rayleigh quotient
        theta = rq(H,K,x);
      end
      u  = U*x;
      Au = AU*x;
      if isempty(B)
        Bu = u;
        uBu = [];
      else
        Bu = BU*x;
        uBu = u'*Bu;
      end
      residual = Au-theta*Bu;
      kmin = 2;
      kmax = min([mindim+1 k]);
    else
      % restart
      found = 0;
      kmin = 1;
      kmax = mindim;
    end
    if nthick == 0
      X = X(:,kmin:kmax);
    elseif found
      X = expand_space(X(:,kmin:kmax-nthick+1), proj([xold(:,2:nthick); zeros(1,nthick-1)], x));
    else
      X = expand_space(X(:,kmin:kmax-nthick), [xold; zeros(1,nthick)]);
    end
    k = size(X,2);
    [U, AU, BU, AtBU] = reduce_dim(X, U, AU, BU, AtBU);
    if found
      AU = projss(AU, Ufound);
      if ~isempty(B)
        BU = projss(BU, Ufound);
      end
      if AtBU_needed
        AtBU = projss(AtBU, Ufound);
      end
    end
    [H, K, UACAU] = reduce_dim_two(X, X, H, K, UACAU);
    if AtBU_needed
      [Q, R, Q1] = reduce_qr(Q, R, X);
      [QAU, QBU] = reduce_dim_two(Q1, X, QAU, QBU);
    end
  end
    
  % Output in Latex format
  if (verbosity == 2) || (verbosity > 2 && ~mod(nr_it, verbosity))
    fprintf('%4d & %4d & %+10.3g %+10.3gi & %.1e & %.2f\n', ...
      nr_it, k, real(theta), imag(theta), normres, relres);
  end
   
  % Solve the correction equation
  if nargout > 4
    nr_prec(nr_it) = nr_prec(nr_it) + ~isempty(M1);
  end

  % If target is not a complex number or residual norm small: take shift = Rayleigh quotient
  small_res = smaller(normres, fix) || ~isnumeric(target);
  if small_res
    shift = theta;
  else
    shift = target;
  end
  
  if isempty(M1)
    if strcmp(expansion, 'jd')           % JD correction equation
      [s, relres, mvs] = ...
        linear_solver(@(x) jd_op(x, A, B, shift, Bu, u, uBu, u, u, [], Ufound, Ufound, Ufound, Ufound, 0), ...
          -residual, solver, maxit_inner, inner_tol);
    elseif strcmp(expansion, 'symjd1')   % Symmetrized JD
%       solver = 'minres';
      solver = 'cr';
      residual = mv(A,residual,1)-theta'*mv(B,residual,1);
      residual = proj(residual, u);
      [s, relres, mvs] = ...
        linear_solver(@(x) jdopti_op(x, A, B, shift, u, Ufound), -residual, solver, maxit_inner, inner_tol);
      mvs = 2*mvs;
    elseif strcmp(expansion, 'symjd2')
      solver = 'lsqr';
      u = mv(A,u,1)-theta'*mv(B,u,1); u = u / norm(u);
      [s, relres, mvs] = ...
        linear_solver(@(x,flag) jdls_op(x, A, B, shift, u, Ufound, flag), -u, solver, maxit_inner, inner_tol);
    elseif strcmp(expansion, 'exact')    % Exact solution
      s = (A-shift*B) \ Bu; % MIQ
      relres = 0;
      mvs = 0;
    elseif strcmp(expansion, 'rqi')      % Accelerated RQI
      [s, relres, mvs] = linear_solver(@(x) mv(A,x)-shift*mv(B,x), Bu, solver, ...
        maxit_inner+1, inner_tol);
      mvs = mvs-1;
    elseif strcmp(expansion, 'resit')    % Residual iteration
      [s, relres, mvs] = linear_solver(@(x) mv(A,x)-shift*mv(B,x), residual, ...
        solver, maxit_inner, inner_tol);
      mvs = mvs-1;
    elseif strcmp(expansion, 'augm')     % Bordered matrix
      [s, relres, mvs] = linear_solver(@(x) bordered(A,x,Bu,shift,u), [-residual; 0], ...
        solver, maxit_inner, inner_tol);
      s = s(1:end-1);
    elseif strcmp(expansion, 'riccati')  % Riccati method
      [V, H1, AV] = krylov(@(x) jd_op(x, A, B, shift, Bu, u, uBu, u, u, [], Ufound, Ufound, Ufound, Ufound, 0), ...
          -residual, maxit_inner, 0);
      [~, ~, s] = extraction_ep(extraction, target, add_row_column([V u], [AV Au], H1));
      s = [V u] * s;
%       [X1, E1] = eig0(add_row_column([V u], [AV Au], H1), [], target);
%       s = [V u] * X1(:,1);
%       s1 = correq_riccati(A, shift, u, -residual, target, Au, maxit_inner);
%       subspace(s,s1)
      mvs = maxit_inner;
    elseif strcmp(expansion, 'krylov')   % Eigenvector from Krylov space
      [V, H1] = krylov(@(x) mv(A,x)-shift*mv(B,x), u, maxit_inner, 0);
      [~, ~, s] = extraction_ep(extraction, target, H1);
      s = V*s;
      mvs = maxit_inner;
    elseif strcmp(expansion, 'expm')     % expm(B\A) operator for rightmost
      [V, H1] = krylov(@(x) minvv(B, mv(A,x)), u, maxit_inner, 0);
      E = eig0(H1, [], 'inf');
      s = V*expm(H1-E(1)*eye(maxit_inner));
      s = s(:,1);
      mvs = maxit_inner;
    else
      error('Unknown expansion method');
    end
    precs = 0;
  else  % Preconditioning
    if strcmp(expansion, 'jd')
      z = prec(M1, M2, Bu);
      Z = [MUfound z];
      [s, relres, mvs, precs] = ...
        linear_solver(@(x) jd_op(x, A, B, shift, Bu, u, uBu, u, u, [], Ufound, Ufound, Ufound, Ufound, 0), ...
          -residual, solver, maxit_inner, inner_tol, ...
             @(x) jd_prec(x, M1, M2, [Ufound u], Z, [Ufound u]'*Z), [], Mtype);
    elseif strcmp(expansion, 'gd')
      s = prec(M1, M2, residual);
      mvs = 0; precs = 1;
    elseif strcmp(expansion, 'gdupdated')
      w = prec(M1, M2, Au);
      w = (w-u)/(u'*w);
      s = jd_prec_updated(residual, M1, M2, u, w, [Ufound u], Z, [Ufound u]'*Z);
      mvs = 0; precs = 1;
    elseif strcmp(expansion, 'gdstab')
      s = prec(M1, M2, theta'*Au+Bu);
      mvs = 0; precs = 1;
    elseif strcmp(expansion, 'gd2')
      s = prec(M1, M2, [Au Bu]);
      mvs = 0; precs = 2;
    elseif strcmp(expansion, 'olsen')     % -M\r + a*M\Bu, where a = u'(M\r) / u'(M\Bu)
      z = prec(M1, M2, Bu);
      Mr = prec(M1, M2, residual);
      s = -Mr + ((u'*Mr) / (u'*z)) * z;
      mvs = 0; precs = 2;
    elseif strcmp(expansion, 'davidson')  % Diagonal preconditioning
      if isempty(B)
        d = diag(A)-shift;
      else
        d = diag(A)-shift*diag(B);
      end
      d(abs(d) < 1e-6) = 1e-6;
      s = -residual ./ d;
      mvs = 0; precs = 0;
    elseif strcmp(expansion, 'subspace')
      [s, relres, mvs, precs] = ...
        linear_solver(@(x) jd_op(x, A, B, shift, U, U, [], U, U, [], Ufound, Ufound, Ufound, Ufound, 0), ...
          -residual, solver, maxit_inner, inner_tol, ...
             @(x) jd_prec(x, M1, M2, [Ufound u], Z, [Ufound u]'*Z), [], Mtype);
    elseif strcmp(expansion, 'updated')
      w = prec(M1, M2, Au);
      w = (w-u)/(u'*w);
      [s, relres, mvs, precs] = ...
        linear_solver(@(x) jd_op(x, A, B, shift, Bu, u, uBu, u, u, [], Ufound, Ufound, Ufound, Ufound, 0), ...
          -residual, solver, maxit_inner, inner_tol, ...
             @(x) jd_prec_updated(x, M1, M2, u, w, [Ufound u], Z, [Ufound u]'*Z), [], Mtype);
    elseif strcmp(expansion, 'jdup')
      [s, relres, mvs, precs] = ...
        linear_solver(@(x) jd_op(x, A, B, shift, Bu, u, uBu, u, u, [], Ufound, Ufound, Ufound, Ufound, 0), ...
          -residual, solver, maxit_inner, inner_tol, ...
             @(x) prec(M1, M2, x), [], Mtype);
    elseif strcmp(expansion, 'jdupm')
      [s, relres, mvs, precs] = ...
        linear_solver(@(x) jd_op(x, A, B, shift, Bu, u, uBu, u, u, [], Ufound, Ufound, Ufound, Ufound, 0), ...
          -residual, solver, maxit_inner, inner_tol, ...
             @(x) proj(prec(M1, M2, x), u), [], Mtype);
    else
      error('Unknown expansion method');
    end
  end
  if nargout > 3
    nr_mv(nr_it)   = nr_mv(nr_it)   + mvs * (1 + ~isempty(B));
  end
  if nargout > 4
    nr_prec(nr_it) = nr_prec(nr_it) + precs*(~isempty(M1));
  end

  % For extra stability (expansion vector s might already be in search space)
  s = projss(s, U);
  if norm(s) < 1e-15
    s = residual; % Expand with residual
    if verbosity > 1
      fprintf('Info: expansion process gives no expansion of search space, taking residual instead\n');
    end
  end
  s = rgs(s, Ufound);
end


if verbosity
  fprintf('**** Exit (%d found) after %d iterations with norm residual = %g\n', length(lambda), nr_it, normres);
end
warning('on', 'MATLAB:rankDeficientMatrix');
