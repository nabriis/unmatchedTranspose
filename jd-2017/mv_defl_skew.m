function v = mv_defl_skew(A, u, U1, V1, U2, V2, transp_flag)

%MV_DEFL_SKEW  Deflated matrix times vector
% function v = mv_defl_skew(A, u, U1, V1, U2, V2, transp_flag)
% In:  V1'*U1 = diag
%      V2'*U2 = diag
% Out: v = (I-U1 V1')A (I-U2 V2')u
%  or: v = (I-V2 U2')A'(I-V1 U1')u  with transpose_flag
%
% See also MV
%
% Revision date: January 22, 2014
% (C) Michiel Hochstenbach 2014

if nargin < 7 || isempty(transp_flag)
  transp_flag = 0;
end

if ~transp_flag
  v = projss_skew(u, U2, V2);
  v = mv(A, v);
  v = projss_skew(v, U1, V1);
else
  v = projss_skew(u, V1, U1);
  v = mv(A, v, transp_flag);
  v = projss_skew(v, V2, U2);
end
