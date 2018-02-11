% This function solves Au=f using SSOR method and returns u,number of
% iterations and 2-norm errors.

function [u,iter,errs] = ssor(A, f, u0, uexact, omega, tol)

f = omega*f;
D = diag(diag(A));
L = tril(A,-1);
R = triu(A,1);
M1 = D + omega*L;               % at n+1/2: D+wL
N1 = M1 - omega*A;              % (1-w)D-wR
M2 = D + omega*R;               % at n+1: D+wR
N2 = M2 - omega*A;              % (1-w)D-wL

errs = zeros(1,1);
iter = 1;
u = u0;
err = norm(u0 - uexact);
errs(1) = err;

while iter < 50000 && err > tol
    u = M1\(N1*u + f);
    u = M2\(N2*u + f);
    err = norm(u - uexact);
    iter = iter + 1;
    errs(iter) = err;
end