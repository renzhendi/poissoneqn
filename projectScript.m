clear

% this script evaluates the solution of the poisson equation given a RHS
% function f and mesh sizes N and M


N=50;
M=50;

%%%%%%%%%%%%% FUNCTION CHOICE %%%%%%%%%%%

RHS = @(x,y) 13*pi^2*sin(2*pi*x).*sin(3*pi*y);
exsolfun = @(x,y) sin(2*pi*x).*sin(3*pi*y);

% RHS = @(x,y) -(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2.*x.^2.*(x-1).^5;
% exsolfun = @(x,y) (x-1).^5.*x.^2.*y.*(y-1);

%%%%%%%%%%%% APPROXIMATE SOLUTION EVALUATION %%%%%%%%%%%%

u=laplaceSolver(N,M,RHS);

error=norm(exsol-u) %display error

%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%

figure()

%-------------approximate solution
subplot(2,1,1)   
surf(u)
title('Approximate Solution')
axis tight

%-------------exact solution
subplot(2,1,2)
exsol=zeros(N+1,M+1);
exsol(2:N,2:M)=exsolfun(X,Y)';
surf(exsol)  %plot exact solution 
title('Exact Solution')
axis tight

%%%%%%%%%%%%%% ERROR %%%%%%%%%%%%%%%%%

error=norm(exsol-u);
disp(error)