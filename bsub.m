% This function implements backward substitution. It solves Ux=y where U is
% an upper triangular matrix, and returns x.

function x = bsub(U,y)

%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%

n = size(U,1);
x = y;
x(n) = y(n)/U(n,n);

%%%%%%%%%%%%%
% main loop %
%%%%%%%%%%%%%

for i = (n-1):-1:1
    x(i) = (y(i) - U(i,i+1:n)*x(i+1:n))/U(i,i);
end