function H = add_row_column(U, V, H)

%ADD_ROW_COLUMN  Add row(s) and column(s) to H = U'*V
% function H = add_row_column(U, V, H)
% In: H: m1 x n1
%     U: .. x m2  (m2 >= m1)
%     V: .. x n2  (n2 >= n1)
%
% Revision date: March 3, 2014
% (C) Michiel Hochstenbach 2014

[m1,n1] = size(H);
m2      = size(U,2);
n2      = size(V,2);
H(m1+1:m2,1:n2) = U(:,m1+1:m2)'*V;
H(1:m2,n1+1:n2) = U'*V(:,n1+1:n2);
% Much faster than H(1:m1,n1+1:n2) = U(:,1:m1)'*V(:,n1+1:n2);
