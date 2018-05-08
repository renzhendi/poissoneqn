% This function solves Ax=b using conjugate gradient algorithm, and returns
% u, number of iterations, and 2-norm errors.
%
% Note: stopping criteria based on 2-norm err, not res.

function [x,iter,errs] = cg_err(A,b,x0,xexact,tol)

errs = zeros(1,10000);
iter = 1;
x = x0;

r = b - A*x0;
err = norm(x0 - xexact);
errs(1) = err;
p = r;

while iter < 10000 && err > tol
    alpha = (p'*r)/(p'*A*p);
    x = x + alpha*p;
    r = b - A*x;
    beta = -(p'*A*r)/(p'*A*p);
    p = r + beta*p;
    err = norm(x - xexact);
    
    iter = iter + 1;
    errs(iter) = err;
end

errs = errs(1:iter);