% This function solves Ax=b using multi-grid method recursively, and returns
% u, number of iterations, and 2-norm errors.
%
% Note: the 2-norm err is computed on the fine grid.

function [u,iter,res_vec] = multigrid(A, b, u0, tol)

u = u0;
iter = 1;
res = norm(b - A*u);
res_vec(iter) = res;

while iter<10000 && res>tol
    u = vcycle(A, b, u);
    res = norm(b - A*u);
    iter = iter+1;
    res_vec(iter) = res;
end