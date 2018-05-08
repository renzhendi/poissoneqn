% This function solves Au=f using SSOR method, and returns u,number of
% iterations, 2-norm errors, and splitted matrices.

function [u,iter,errs,M,N] = ssor_err(A, f, u0, uexact, omega, tol)

f = omega*f;
D = diag(diag(A));
L = tril(A,-1);
R = triu(A,1);
M1 = D + omega*L;               % at n+1/2: D+wL
N1 = M1 - omega*A;              % (1-w)D-wR
M2 = D + omega*R;               % at n+1: D+wR
N2 = M2 - omega*A;              % (1-w)D-wL
M = M1*M2/(omega*(2-omega)*A(1,1));
N = M - A;

errs = zeros(1,10000);
iter = 1;
u = u0;
err = norm(u0 - uexact);
errs(1) = err;

while iter < 10000 && err > tol
    u = M1\(N1*u + f);
    u = M2\(N2*u + f);
    err = norm(u - uexact);
    iter = iter + 1;
    errs(iter) = err;
end

errs = errs(1:iter);