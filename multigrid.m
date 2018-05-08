% This function solves Ax=b using multi-grid method recursively, and returns
% u, number of iterations (v-cycles), and 2-norm residuals.
%
% Note: the 2-norm res is computed on the fine grid.

function [u,iter,res_vec] = multigrid(A, b, u0, tol)

res_vec = zeros(1,10000);
iter = 1;
u = u0;
res = norm(b - A*u);
res_vec(1) = res;

while iter < 10000 && res>tol
    u = vcycle(A, b, u);
    res = norm(b - A*u);
    iter = iter+1;
    res_vec(iter) = res;
end

res_vec = res_vec(1:iter);