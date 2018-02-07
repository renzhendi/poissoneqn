clear
clc
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script solves general Poisson's equations numerically. It allows users
% to specify M and N (size of the grid), and input f (the RHS of Poisson's
% equation), Then Au=f is solved using the simple finite difference scheme.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
solverIndex = 1;
timingBoolean = 0;
initGuessType = 2;
relaxation = 1;
savePlotBoolean = timingBoolean;

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
N = 100;                  % number of intervals in x direction
M = 50;                  % number of intervals in y direction
A_name = sprintf('matA_m%d_n%d.mat',M,N);
if exist(A_name,'file') == 2
    load(A_name);
else
    A = createA(M,N);
end

%%%%%%%%%%%%%%%%%%%%%%
% RHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
x = (1:N-1)/N;                          % horizontal axis
y = ((1:M-1)/M)';                       % vertical axis
f1=@(x,y)13*pi^2*sin(2*pi*x).*sin(3*pi*y);
f1mat = f1(x,y)';                       % (N-1)*(M-1) matrix
f1vec = f1mat(:);                       % (N-1)(M-1) vector
f2=@(x,y)-(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5;
f2mat = f2(x,y)';                       % (N-1)*(M-1) matrix
f2vec = f2mat(:);                       % (N-1)(M-1) vector

%%%%%%%%%%%%%%%%%%
% numerical soln %
%%%%%%%%%%%%%%%%%%
[u1,t1,iter1,errs1] = solvers(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation);   % (M-1)(N-1) vector and avg run time
u1mat = vec2mat(u1,N-1);                                   % (M-1)*(N-1) matrix
[u2,t2,iter2,errs2] = solvers(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation);
u2mat = vec2mat(u2,N-1);

%%%%%%%%%%%%%%
% exact soln %
%%%%%%%%%%%%%%
s1=@(x,y)sin(2*pi*x).*sin(3*pi*y);
s1mat = s1(x,y);
s2=@(x,y)(x-1).^5.*x.^2.*y.*(y-1);
s2mat = s2(x,y);

%%%%%%%%%%%%%%%%
% verification %
%%%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1)    % plot the numerical solution
s11 = surf(x,y,u1mat,'FaceAlpha',0.7); s11.EdgeColor = 'none';
title(sprintf('Eqn1: M=%d, N=%d completed in %0.6f secs',M,N,t1));
axis([0 1 0 1 -1 1]); colorbar;
xlabel('x'); ylabel('y'); zlabel('numerical soln');
subplot(3,1,2)    % plot the exact solution
s12 = surf(x,y,s1mat,'FaceAlpha',0.7); s12.EdgeColor = 'none';
axis([0 1 0 1 -1 1]); colorbar;
xlabel('x'); ylabel('y'); zlabel('exact soln');
subplot(3,1,3)    % plot the residual (numerical - exact)
s13 = surf(x,y,u1mat-s1mat,'FaceAlpha',0.7); s13.EdgeColor = 'none';
colorbar; xlabel('x'); ylabel('y'); zlabel('numerical - exact');
if savePlotBoolean
    print(sprintf('plots/Eqn1_bs_m%d_n%d.png', M, N),'-dpng');
end

figure(2)
subplot(3,1,1)
s21 = surf(x,y,u2mat,'FaceAlpha',0.7); s21.EdgeColor = 'none';
title(sprintf('Eqn2: M=%d, N=%d completed in %0.6f secs',M,N,t2));
axis([0 1 0 1 -0.001 0.004]); colorbar;
xlabel('x'); ylabel('y'); zlabel('numerical soln');
xlabel('x'); ylabel('y'); zlabel('numerical soln');
subplot(3,1,2)
s22 = surf(x,y,s2mat,'FaceAlpha',0.7); s22.EdgeColor = 'none';
axis([0 1 0 1 -0.001 0.004]); colorbar;
xlabel('x'); ylabel('y'); zlabel('exact soln');
subplot(3,1,3)
s23 = surf(x,y,u2mat-s2mat,'FaceAlpha',0.7); s23.EdgeColor = 'none';
colorbar; xlabel('x'); ylabel('y'); zlabel('numerical - exact');
if savePlotBoolean
    print(sprintf('plots/Eqn2_bs_m%d_n%d.png', M, N),'-dpng');
end