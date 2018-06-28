function [w] = est_step_size_mf(A,B,alpha2,factor)
%EST_STEP_SIZE estimate of step size using full matrix

%Get minimum real part of eigenvalues of BA as regularization parameter
%A_m = A*speye(size(A,2));
%B_m = B*speye(size(B,2));
opts.maxit = 500; opts.tol = 10^(-2);
eigenvalues_full = eigs(A*B,size(A,1),'lr',opts);

%Find the smallest bound to select step size (w)

for i = 1:length(eigenvalues_full)
    
    numerator = alpha2+real(eigenvalues_full(i));
    denominator = real(eigenvalues_full(i))^2+imag(eigenvalues_full(i))^2+alpha2*(alpha2+2*real(eigenvalues_full(i)));
    
    bound(i) = factor * (numerator/denominator); 
    
end

[w,idx] = min(bound);

end

