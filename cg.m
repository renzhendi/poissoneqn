% This function solves Ax=b using conjugate gradient algorithm, and returns
% u, number of iterations, and 2-norm errors.
%
% Note: stopping criteria based on 2-norm err, not res.

function [x,iter,res_vec] = cg(A, b, x0, tol)

res_vec = zeros(1,50000);
iter = 1;
x = x0;

r = b - A*x0;
res = norm(r);
res_vec(1) = res;
p = r;

while iter < 50000 && res > tol
    alpha = (p'*r)/(p'*A*p);
    x = x + alpha*p;
    r = b - A*x;
    beta = -(p'*A*r)/(p'*A*p);
    p = r + beta*p;
    res = norm(r);
    
    iter = iter + 1;
    res_vec(iter) = res;
end

res_vec = res_vec(1:iter);