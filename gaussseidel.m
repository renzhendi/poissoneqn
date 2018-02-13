% This function solves Au=f using Gauss-Seidel/SOR method and returns u,
% number of iterations and 2-norm errors.

function [u,iter,errs,M,N] = gaussseidel(A, f, u0, uexact, omega, tol)

if omega == 1
    M = tril(A);                          % M = D+L
    N = M - A;                            % N = -R
else
    f = omega*f;
    M = diag(diag(A)) + omega*tril(A,-1); % M = D+wL
    N = M - omega*A;                      % N = (1-w)D-wR
end

errs = zeros(1,1);
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