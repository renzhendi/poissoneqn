% This function solves Au=f using (regular/relaxed) Jacobi's method returns
% u, number of iterations and 2-norm errors.
%
% Component-wise implementation.

function [u,iter,errs] = jacobi2(A, f, u0, uexact, theta, tol)

M = diag(diag(A));                   % M = D
d = M(1,1);

if theta == 1
    N = M - A;                       % N = -(L+R)
else
    f = theta*f;
    N = (1-theta)*M+theta*(M-A);     % N = (1-t)D-t(L+R)
end

errs = zeros(1,1);
iter = 1;
u = u0;
err = norm(u0 - uexact);
errs(1) = err;

while iter < 50000 && err > tol
    u = (N*u + f)/d;
    err = norm(u - uexact);
    iter = iter + 1;
    errs(iter) = err;
end