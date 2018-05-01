% This function solves Ax=b using conjugate gradient algorithm, and returns
% u, number of iterations, 2-norm errors, A-norm errors and 2-norm residuals.
%
% Note: stopping criteria based on 2-norm err, not res.

function [x,iter,errs,errsA,res] = cg_reserr(A, b, x0, xexact, tol)

res = zeros(1,1);
errsA = res;
errs = res;
iter = 1;
x = x0;

r = b - A*x0;
res(1) = norm(r);
e = x0 - xexact;
errsA(1) = sqrt(e'*A*e);
err = norm(e);
errs(1) = err;
p = r;

while iter < 50000 && err > tol
    alpha = (p'*r)/(p'*A*p);
    x = x + alpha*p;
    r = b - A*x;
    beta = -(p'*A*r)/(p'*A*p);
    p = r + beta*p;
    e = x - xexact;
    
    iter = iter + 1;
    res(iter) = norm(r);
    errsA(iter) = sqrt(e'*A*e);
    err = norm(e);
    errs(iter) = err;
end