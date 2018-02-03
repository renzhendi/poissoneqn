N=100;
M=100;

hx=1/N;
hy=1/M;
A=Amatrix2(N,M);
x=(1:N-1)*hx;
y=(1:M-1)*hy;

[X,Y]=meshgrid(x,y);

%f=13*pi^2.*sin(2*pi*X).*sin(3*pi*Y);
f=13*pi^2*sin(2*pi*x)'*sin(3*pi*y);

fvec=reshape(f,(N-1)*(M-1),1);

uvec=A\fvec;
u=zeros(N+1,M+1);
u(2:N,2:M)=reshape(uvec,N-1,M-1);

%plot approximate solution
figure(1)
surf(u)
title('Approximate Solution')
axis tight

%plot exact solution 
figure(2)
exsol=zeros(N+1,M+1);
exsol(2:N,2:M)=sin(2*pi*x)'*sin(3*pi*y);
surf(exsol)
title('Exact Solution')
axis tight

error=norm(exsol-u)