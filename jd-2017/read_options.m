function varargout = read_options(opts, varargin)

%READ_OPTIONS  Read options from a structure
% function varargout = read_options(opts, varargin)
%
% Most frequently used fields of option:                             Default
%   m, n         = size(A) when A is a function                      size(A)
%   nr           number of desired pairs or triples                  1
%   target       complex number | '-inf' | 'inf' | 'x-axis' | 'abs'  'abs'
%   u1           initial left  space (may be more-dim.)              randn(m,1)
%   v1           initial right space (dim(v1)=dim(u1))               randn(n,1)
%   x0           starting vector for Ax=b                            []
%   rand(n,1)
%   tol          tolerance of the outer iteration                    1e-6
%   absrel       absolute or relative tolerance                      'rel'
%   method       choice of method for certain algorithms             1
%   mindim       minimum dimension of subspaces                      30
%   maxdim       maximum dimension of subspaces                      60
%   maxit        maximum number of outer iterations                  5000
%   expansion    method to compute expansions                        'jd'
%   solver       solver for inner iterations                         'gmres'
%   maxit_inner  maximum number of inner iterations                  10
%   inner_tol    tolerance of the inner iteration                    0
%   expansion    how to expand the search space                      'arnoldi'
%   extraction   how to perform extraction from a search space       'standard'
%   fix          fix target until residual norm < fix                0.01
%   M1, M2       M=M1*M2: preconditioner in factored form            [], []
%   Mtype        left | right | implicit preconditioning             'left'
%   thick        thick restart                                       1
%   krylov       start method with Krylov space of dimension maxdim  1
%   block        block size                                          1
%   verbosity    output desired (0..2)                               0
% Less frequently used fields:
%   L1, L2       L=L1*L2: preconditioner in factored form (for MEP)  [], []
%   u2, v2       initial spaces for 2nd equation (for MEP etc)       [], []
%   full         full orthogonalization in (two-sided) Lanczos       1
%   corrected    corrected approximation for matrix functions        0
%   threshold    switch to harmonic extraction if ||r|| <            1
%   value_extr   from vector for polynomial eigenvalue problem       []
%
% Revision date: March 18, 2014
% (C) Michiel Hochstenbach 2014

if isempty(opts)
  names = [];
else
  names = fieldnames(opts);
end
for j = 1:nargin-1
  argin = varargin{j};
  I = strmatch(argin, names, 'exact');
  if ~isempty(I)
%     field = getfield(opts, names{I});
    field = opts.(names{I});
  else
    field = [];
  end
%   if isfield(opts, argin)
%     field = getfield(opts, argin);
%   else
%     field = [];
%   end
  if isempty(field) 
    switch(argin)         % assign default value
      case 'nr',          field = 1;
      case 'u1',          field = qr_pos(randn(opts.m, getfield0(opts, 'block', 1)),0);
      case 'v1',          field = qr_pos(randn(opts.n, getfield0(opts, 'block', 1)),0);
      case 'x0'
      case 'tol',         field = 1e-6;
      case 'absrel',      field = 'rel';
      case 'method',      field = 1;
      case 'mindim',      field = 30;
      case 'maxdim',      field = 60;
      case 'maxit',       field = 5000;
      case 'expansion',   field = 'jd';
      case 'solver',      field = 'gmres';
      case 'maxit_inner', field = 10;
      case 'inner_tol',   field = 0;
      case 'target',      field = 'abs';
      case 'extraction',  field = 'standard';
      case 'fix',         field = 1e-2;
      case 'M1'
      case 'M2'
      case 'Mtype',       field = 'left';
      case 'krylov',      field = 1;
      case 'thick',       field = 1;
      case 'block',       field = 1;
      case 'verbosity',   field = 0;
      case 'L1',
      case 'L2',
      case 'corrected',   field = 0;
      case 'threshold',   field = 1;
      case 'value_extr',  field = [];
      case 'full',        field = 1;
      case 'u2',          field = randn1(opts.n2);
      case 'v2',          field = randn1(opts.n2);
      otherwise,          field = [];
    end
  end
  varargout(j) = {field};
end
