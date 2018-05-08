% This function solves Au=f using SSOR method, and returns u,number of
% iterations, 2-norm errors, and splitted matrices.

function [u,iter,res_vec,M,N] = ssor(A, f, u0, omega, tol)

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

res_vec = zeros(1,10000);
iter = 1;
u = u0;
res = norm(f - A*u);
res_vec(1) = res;

while iter < 10000 && res > tol
    u = M1\(N1*u + f);
    u = M2\(N2*u + f);
    res = norm(f - A*u);
    iter = iter + 1;
    res_vec(iter) = res;
end

res_vec = res_vec(1:iter);