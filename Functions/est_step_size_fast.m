function [w] = est_step_size_fast(A,B,alpha2,factor)
%EST_STEP_SIZE_FAST lower bound estimate of step size

opts.maxit = 1000; opts.tol = 10^(-2);
eigenvalues_sr = eigs(A*B,1,'sr',opts);
eigenvalues_lr = eigs(A*B,1,'lr',opts);
eigenvalues_lm = eigs(A*B,1,'lm',opts);

w = factor*(eigenvalues_sr+alpha2)/(eigenvalues_lm+alpha2*(alpha2+2*eigenvalues_lr));

end

