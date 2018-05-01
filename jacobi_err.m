% This function solves Au=f using (regular/relaxed) Jacobi's method, and
% returns u, number of iterations, 2-norm errors, and splitted matrices.

function [u,iter,errs,M,N] = jacobi_err(A, f, u0, uexact, theta, tol)

M = diag(diag(A));                   % M = D

if theta == 1
    N = M - A;                       % N = -(L+R)
else
    f = theta*f;
    N = (1-theta)*M+theta*(M-A);     % N = (1-t)D-t(L+R)
end

errs = zeros(1,50000);
iter = 1;
u = u0;
err = norm(u0 - uexact);
errs(1) = err;

while iter < 50000 && err > tol
    u = M\(N*u + f);
    err = norm(u - uexact);
    iter = iter + 1;
    errs(iter) = err;
end

errs = errs(1:iter);