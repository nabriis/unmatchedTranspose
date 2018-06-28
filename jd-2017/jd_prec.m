function y = jd_prec(x, M1, M2, Mu, v, vMu)

%JD_PREC  Perform action of JD preconditioning in correction equation
% function y = jd_prec(x, M1, M2, Mu, v, vMu)
%
% The inverse of (I-u (w'u)\w')M: v_perp -> w_perp is
%                (I-M^{-1}u (v'M^{-1}u)\v') M^{-1}
% Mu = M^{-1}u,  vMu = v'M^{-1}u
% y = (I-Mu (v'Mu)\v') M^{-1}x
%
% See also JD, JDSYM, JDPEP, PREC, PROJSS
%
% Revision date: June 3, 2016
% (C) Michiel Hochstenbach 2016

y = projss(prec(M1, M2, x), Mu, v, vMu);
