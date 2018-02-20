% This function solves Ax=b using multi-grid method recursively, and returns
% u, number of iterations, and 2-norm errors.
%
% Note: the 2-norm err is computed on the fine grid.

function [x,iter,errs] = multigrid(A, b, x0, xexact, tol)

M = diag(diag(A));                   % M = D
N = M - A;                       % N = -(L+R)

errs = zeros(1,1);
iter = 1;
x = x0;

I = speye(size(A,1)/2);
R = 0.5*kron(I,[1 1]);
P = 2*R';
Ac = R*A*P;

r = b - A*x0;
err = norm(x0 - xexact);
errs(1) = err;

while iter < 50000 && err > tol
    rc = R*r;
    ec = Ac\rc;     % should call multigrid recursively
    e = P*ec;
    x = x+e;
    r = b - A*x;
    err = norm(x - xexact);
    
    iter = iter + 1;
    errs(iter) = err;
end