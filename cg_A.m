% This function solves Ax=b using conjugate gradient algorithm, and returns
% u, number of iterations, 2-norm errors, A-norm errors and 2-norm residuals.
%
% Note: stopping criteria based on 2-norm err, not res.

function [x,iter,errA,res_vec] = cg_A(A,b,x0,tol)

res_vec = zeros(1,50000);
errA = res_vec;
iter = 1;
xexact = A\b;
x = x0;

r = b - A*x0;
res = norm(r);
res_vec(1) = res;
e = x0 - xexact;
errA(1) = sqrt(e'*A*e);
p = r;

while iter < 50000 && res > tol
    alpha = (p'*r)/(p'*A*p);
    x = x + alpha*p;
    r = b - A*x;
    beta = -(p'*A*r)/(p'*A*p);
    p = r + beta*p;
    e = x - xexact;
    res = norm(r);
    
    iter = iter + 1;
    errA(iter) = sqrt(e'*A*e);
    res_vec(iter) = res;
end