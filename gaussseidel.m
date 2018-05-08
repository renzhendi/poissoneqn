% This function solves Au=f using Gauss-Seidel/SOR method, and returns u,
% number of iterations, 2-norm errors, and splitted matrices.

function [u,iter,res_vec,M,N] = gaussseidel(A, f, u0, omega, tol)

if omega == 1
    M = tril(A);                          % M = D+L
    N = M - A;                            % N = -R
else
    f = omega*f;
    M = diag(diag(A)) + omega*tril(A,-1); % M = D+wL
    N = M - omega*A;                      % N = (1-w)D-wR
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