function [u,exsol] = laplaceSolver(N,M,RHS,exsolfun)
% function that solves poisson equation using multigrid method and 
% backslash given a RHS function

% inputs: N: length of vector x
%         M: length of vector y
%         RHS: function f such that Au=f

hx=1/N;
hy=1/M;

%construct matrix A
C=-1/hy^2*eye(N-1,N-1);
B=toeplitz([2/hx^2+2/hy^2,-1/hx^2,zeros(1,N-3)]);
A=kron(eye(M-1),B)+kron(toeplitz([0,1,zeros(1,M-3)]),C);

%construct vectors x,y
x=(1:N-1)*hx;
y=(1:M-1)*hy;

[X,Y]=meshgrid(x,y);

exsol=zeros(N+1,M+1);
exsol(2:N,2:M)=exsolfun(X,Y)';

f=RHS(X,Y)';
fvec=reshape(f,(N-1)*(M-1),1); 

%solve for u 
uvec=A\fvec;
u=zeros(N+1,M+1);
u(2:N,2:M)=reshape(uvec,N-1,M-1);



