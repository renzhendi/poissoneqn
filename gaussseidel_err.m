% This function solves Au=f using Gauss-Seidel/SOR method, and returns u,
% number of iterations, 2-norm errors, and splitted matrices.

function [u,iter,errs,M,N] = gaussseidel_err(A, f, u0, uexact, omega, tol)

if omega == 1
    M = tril(A);                          % M = D+L
    N = M - A;                            % N = -R
else
    f = omega*f;
    M = diag(diag(A)) + omega*tril(A,-1); % M = D+wL
    N = M - omega*A;                      % N = (1-w)D-wR
end

errs = zeros(1,10000);
iter = 1;
u = u0;
err = norm(u0 - uexact);
errs(1) = err;

while iter < 10000 && err > tol
    u = M\(N*u + f);
    err = norm(u - uexact);
    iter = iter + 1;
    errs(iter) = err;
end

errs = errs(1:iter);