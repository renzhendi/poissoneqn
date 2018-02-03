N=100;
h=1/N;
A=Amatrix(N);
x=(1:N-1)*h;
y=(1:N-1)*h;

f=13*pi^2*sin(2*pi*x)'*sin(3*pi*y);
fvec=reshape(f,(N-1)^2,1);

uvec=A\fvec;
u=zeros(N+1,N+1);
u(2:N,2:N)=reshape(uvec,N-1,N-1);

%plot approximate solution
figure(1)
surf(u)
title('Approximate Solution')
axis tight

%plot exact solution 
figure(2)
exsol=zeros(N+1,N+1);
exsol(2:N,2:N)=sin(2*pi*x)'*sin(3*pi*y);
surf(exsol)
title('Exact Solution')
axis tight

error=norm(exsol-u);