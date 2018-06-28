function [x, flag, relres, iter, resvec] = gmres0(A, b, tol, maxit, M1, M2, x)

%GMRES   Generalized Minimum Residual Method
% function [x, flag, relres, iter, resvec] = gmres0(A, b, tol, maxit, M1, M2, x)
%
% $Revision: 1.21.4.15 $ $Date: 2011/05/17 02:33:07 $
% Adapted by Michiel Hochstenbach, March 20, 2014

% Assign default values to unspecified parameters
if nargin < 3 || isempty(tol)
  tol = 1e-6;
end
if nargin < 4 || isempty(maxit)
  maxit = 10;
end
if nargin < 5
  M1 = [];
end
if nargin < 6
  M2 = [];
end
if nargin < 7
  x = [];
end
n = length(b);

n2b = sqrt(b'*b);                  % Norm of rhs vector, b

% Set up for the method
flag = 1;
xmin = x;                        % Iterate which has minimal residual so far
imin = 0;                        % "Outer" iteration at which xmin was computed
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolb = tol * n2b;                % Relative tolerance
evalxm = 0;
stag = 0;
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;
minupdated = 0;

if isempty(x)
  r = b;
else
  r = b - mv(A,x);
end
normr = sqrt(r'*r);              % Norm of initial residual
minv_b = b;

if ~isempty(M1)
  r = minvv(M1,r);
  if ~isempty(x)
    minv_b = minvv(M1,b);
  else
    minv_b = r;
  end
end

if ~isempty(M2)
  r = minvv(M2,r);
  if ~isempty(x)
    minv_b = minvv(M2,minv_b);
  else
    minv_b = r;
  end
end

normr = sqrt(r'*r);               % norm of the preconditioned residual
n2minv_b = sqrt(minv_b'*minv_b);  % norm of the preconditioned rhs
clear minv_b;
tolb = tol * n2minv_b;

resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = normr;                % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin

%  Preallocate J to hold the Given's rotation constants.
J = zeros(2,inner);

U = zeros(n,inner);
R = zeros(inner,inner);
w = zeros(inner+1,1);

for outiter = 1 : outer
  %  Construct u for Householder reflector.
  %  u = r + sign(r(1))*||r||*e1
  u = r;
  normr = sqrt(r'*r);
  beta = scalarsign(r(1))*normr;
  u(1) = u(1) + beta;
  u = u * (1 / sqrt(u'*u));

  U(:,1) = u;

  %  Apply Householder projection to r.
  %  w = r - 2*u*u'*r;
  w(1) = -beta;

  for initer = 1 : inner
    %  Form P1*P2*P3...Pj*ej.
    %  v = Pj*ej = ej - 2*u*u'*ej
    v = -2*(u(initer)')*u;
    v(initer) = v(initer) + 1;
    %  v = P1*P2*...Pjm1*(Pj*ej)
    for k = (initer-1):-1:1
      v = v - U(:,k)*(2*(U(:,k)'*v));
    end
    %  Explicitly normalize v to reduce the effects of round-off.
    v = v * (1 / sqrt(v'*v));

    %  Apply A to v.
    v = mv(A,v);
    %  Apply Preconditioner.
    if ~isempty(M1)
      v = minvv(M1,v);
    end

    if ~isempty(M2)
      v = minvv(M2,v);
    end
    %  Form Pj*Pj-1*...P1*Av.
    for k = 1:initer
      v = v - U(:,k)*(2*(U(:,k)'*v));
    end

    %  Determine Pj+1.
    if initer ~= length(v)
      %  Construct u for Householder reflector Pj+1.
      u = [zeros(initer,1); v(initer+1:end)];
      alpha = sqrt(u'*u);
      if alpha ~= 0
        alpha = scalarsign(v(initer+1))*alpha;
        %  u = v(initer+1:end) +
        %        sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
        u(initer+1) = u(initer+1) + alpha;
        u = u * (1 / sqrt(u'*u));
        U(:,initer+1) = u;

        %  Apply Pj+1 to v.
        %  v = v - 2*u*(u'*v);
        v(initer+2:end) = 0;
        v(initer+1) = -alpha;
      end
    end

    %  Apply Given's rotations to the newly formed v.
    for colJ = 1:initer-1
      tmpv = v(colJ);
      v(colJ)   = conj(J(1,colJ))*v(colJ) + conj(J(2,colJ))*v(colJ+1);
      v(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*v(colJ+1);
    end

    %  Compute Given's rotation Jm.
    if ~(initer==length(v))
      rho = norm(v(initer:initer+1));
      J(:,initer) = v(initer:initer+1)./rho;
      w(initer+1) = -J(2,initer).*w(initer);
      w(initer) = conj(J(1,initer)).*w(initer);
      v(initer) = rho;
      v(initer+1) = 0;
    end

    R(:,initer) = v(1:inner);

    normr = abs(w(initer+1));
    resvec((outiter-1)*inner+initer+1) = normr;
    normr_act = normr;

    if normr <= tolb || stag >= maxstagsteps || moresteps
      if evalxm == 0
        ytmp = R(1:initer,1:initer) \ w(1:initer);
        additive = U(:,initer)*(-2*ytmp(initer)*conj(U(initer,initer)));
        additive(initer) = additive(initer) + ytmp(initer);
        for k = initer-1 : -1 : 1
          additive(k) = additive(k) + ytmp(k);
          additive = additive - U(:,k)*(2*(U(:,k)'*additive));
        end
        if sqrt(additive'*additive) < eps*sqrt(x'*x)
          stag = stag + 1;
        else
          stag = 0;
        end
        xm = x + additive;
        evalxm = 1;
      elseif evalxm == 1
        addvc = [-(R(1:initer-1,1:initer-1)\R(1:initer-1,initer))*...
          (w(initer)/R(initer,initer)); w(initer)/R(initer,initer)];
        if sqrt(addvc'*addvc) < eps*sqrt(xm'*xm)
          stag = stag + 1;
        else
          stag = 0;
        end
        additive = U(:,initer)*(-2*addvc(initer)*conj(U(initer,initer)));
        additive(initer) = additive(initer) + addvc(initer);
        for k = initer-1 : -1 : 1
          additive(k) = additive(k) + addvc(k);
          additive = additive - U(:,k)*(2*(U(:,k)'*additive));
        end
        xm = xm + additive;
      end
      r = b - mv(A,xm);
      if sqrt(r'*r) <= tol*n2b
        x = xm;
        flag = 0;
        iter = [outiter, initer];
        break
      end
      minv_r = r;
      if ~isempty(M1)
        minv_r = minvv(M1,r);
      end
      if ~isempty(M2)
        minv_r = minvv(M2,minv_r);
      end

      normr_act = sqrt(minv_r'*minv_r);
      resvec((outiter-1)*inner+initer+1) = normr_act;

      if normr_act <= normrmin
        normrmin = normr_act;
        imin = outiter;
        jmin = initer;
        xmin = xm;
        minupdated = 1;
      end

      if normr_act <= tolb
        x = xm;
        flag = 0;
        iter = [outiter, initer];
        break
      else
        if stag >= maxstagsteps && moresteps == 0
          stag = 0;
        end
        moresteps = moresteps + 1;
        if moresteps >= maxmsteps
          if ~warned
            warning(message('MATLAB:gmres:tooSmallTolerance'));
          end
          flag = 3;
          iter = [outiter, initer];
          break;
        end
      end
    end

    if normr_act <= normrmin
      normrmin = normr_act;
      imin = outiter;
      jmin = initer;
      minupdated = 1;
    end

    if stag >= maxstagsteps
      flag = 3;
      break;
    end
  end         % ends inner loop

  evalxm = 0;

  if flag ~= 0
    if minupdated
      idx = jmin;
    else
      idx = initer;
    end
    y = R(1:idx,1:idx) \ w(1:idx);
    additive = U(:,idx)*(-2*y(idx)*conj(U(idx,idx)));
    additive(idx) = additive(idx) + y(idx);
    for k = idx-1 : -1 : 1
      additive(k) = additive(k) + y(k);
      additive = additive - U(:,k)*(2*(U(:,k)'*additive));
    end
    x = x + additive;
    xmin = x;
    r = b - mv(A,x);
    minv_r = r;
    if ~isempty(M1)
      minv_r = minvv(M1,r);
      if ~all(isfinite(minv_r))
        flag = 2;
        break
      end
    end
    if ~isempty(M2)
      minv_r = minvv(M2,minv_r);
      if ~all(isfinite(minv_r))
        flag = 2;
        break
      end
    end
    normr_act = sqrt(minv_r'*minv_r);
    r = minv_r;
  end

  if normr_act <= normrmin
    xmin = x;
    normrmin = normr_act;
    imin = outiter;
    jmin = initer;
  end

  if flag == 3
    break;
  end
  if normr_act <= tolb
    flag = 0;
    iter = [outiter, initer];
    break;
  end
  minupdated = 0;
end         % ends outer loop

% Returned solution is that with minimum residual
if flag == 0
  relres = normr_act / n2minv_b;
else
  x = xmin;
  iter = [imin jmin];
  relres = normr_act / n2minv_b;
end

% Truncate the zeros from resvec
if nargout > 4
  if flag <= 1 || flag == 3
    resvec = resvec(1:(outiter-1)*inner+initer+1);
    indices = resvec==0;
    resvec = resvec(~indices);
  else
    if initer == 0
      resvec = resvec(1:(outiter-1)*inner+1);
    else
      resvec = resvec(1:(outiter-1)*inner+initer);
    end
  end
end

function sgn = scalarsign(d)
sgn = sign(d);
if sgn == 0
  sgn = 1;
end
