% This function solves Au=f using (regular/relaxed) Jacobi's method, and
% returns u, number of iterations, 2-norm errors, and splitted matrices.

function [u,iter,res_vec,M,N] = jacobi(A, f, u0, theta, tol)

M = diag(diag(A));                   % M = D

if theta == 1
    N = M - A;                       % N = -(L+R)
else
    f = theta*f;
    N = (1-theta)*M+theta*(M-A);     % N = (1-t)D-t(L+R)
end

res_vec = zeros(1,10000);
iter = 1;
u = u0;
res = norm(f - A*u);
res_vec(1) = res;

while iter < 10000 && res > tol
    u = M\(N*u + f);
    res = norm(f - A*u);
    iter = iter + 1;
    res_vec(iter) = res;
end

res_vec = res_vec(1:iter);