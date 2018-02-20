% This function implements forward substitution. It solves Ly=b where L is
% a lower triangular matrix, and returns y.

function y = fsub(L,b)

%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%

n = size(L,1);
y = b;
y(1) = b(1)/L(1,1);

%%%%%%%%%%%%%
% main loop %
%%%%%%%%%%%%%

for i = 2:n
    y(i) = (b(i) - L(i,1:i-1)*y(1:i-1))/L(i,i);
end