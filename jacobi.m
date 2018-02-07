% This function solves Au=f using (regular/relaxed) Jacobi's method and
% returns u and number of iterations.

function [u,iter,errs] = jacobi(A, f, u0, uexact, theta)

M = diag(diag(A));               % M = D
N = M - A;                       % N = -(L+R)

if theta ~= 1
    f = theta*f;
    N = (1-theta)*M+theta*N;     % N = (1-t)D-t(L+R)
end

errs = zeros(10000,1);
iter = 1;
u = u0;
err = norm(u0 - uexact);
errs(1) = err;

while iter < 20000 && err > 10^(-8)
    u = (N*u + f)./diag(M);
    err = norm(u - uexact);
    iter = iter + 1;
    errs(iter) = err;
end

errs = errs(1:iter);