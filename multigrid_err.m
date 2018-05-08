% This function solves Ax=b using multi-grid method recursively, and returns
% u, number of iterations (v-cycles), and 2-norm errors.
%
% Note: the 2-norm err is computed on the fine grid.

function [u,iter,errs] = multigrid_err(A,b,u0,uexact,tol)

errs = zeros(1,10000);
iter = 1;
u = u0;

err = norm(u - uexact);
errs(iter) = err;

while iter < 10000 && err>tol
    u = vcycle(A, b, u);
    err = norm(u - uexact);
    iter = iter+1;
    errs(iter) = err;
end

errs = errs(1:iter);