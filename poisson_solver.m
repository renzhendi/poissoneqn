clear
clc
% This script solves general Poisson's equations numerically. It allows
% users to input f, the RHS of Poisson's equation, and solves Au=f using
% the simple finite difference scheme.

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
N = 1000;                                                 % x, number of intervals in discretization
M = 500;                                                 % y, number of intervals
h1 = 1/N;                                                % x, interval size
h2 = 1/M;                                                % y, interval size
diag0 = ones(N-1,1)*2*(h1^2+h2^2)/(h1^2*h2^2);           % trivial step in constructing B
diag1 = -ones(N-1,1)/(h1^2);
B = spdiags([diag1, diag0, diag1], -1:1, N-1, N-1);      % block B, (N-1)*(N-1) matrix
C = -speye(N-1)/(h2^2);                                  % block C, (N-1)*(N-1) matrix
mat_cells = {B,C,zeros(N-1)};                            % trivial step in constructing A
A = cell2mat(mat_cells(toeplitz([1,2,ones(1,M-3)*3])));  % mat A, (M-1)(N-1)*(M-1)(N-1) matrix

%%%%%%%%%%%%%%%%%%%%%%
% RHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
x = (1:N-1)*h1;                          % horizontal axis
y = ((1:M-1)*h2)';                       % vertical axis
f1=@(x,y)13*pi^2*sin(2*pi*x).*sin(3*pi*y);
f1mat = f1(x,y)';                        % (N-1)*(M-1) matrix
f1vec = f1mat(:);                        % (N-1)(M-1) vector
f2=@(x,y)-(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5;
f2mat = f2(x,y)';                        % (N-1)*(M-1) matrix
f2vec = f2mat(:);                        % (N-1)(M-1) vector

%%%%%%%%%%%%%%%%%%
% numerical soln %
%%%%%%%%%%%%%%%%%%
u1 = A\f1vec;                % backslash solver (Gaussian)
u1mat = vec2mat(u1,N-1);     % (M-1)*(N-1) matrix
u2 = A\f2vec;
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
axis([0 1 0 1 -1 1]); title('Numerical Solution'); colorbar;
subplot(3,1,2)    % plot the exact solution
s12 = surf(x,y,s1mat,'FaceAlpha',0.7); s12.EdgeColor = 'none';
axis([0 1 0 1 -1 1]); title('Exact Solution'); colorbar
subplot(3,1,3)    % plot the residual (numerical - exact)
s13 = surf(x,y,u1mat-s1mat,'FaceAlpha',0.7); s13.EdgeColor = 'none';
title('Numerical - Exact'); colorbar;

figure(2)
subplot(3,1,1)
s21 = surf(x,y,u2mat,'FaceAlpha',0.7); s21.EdgeColor = 'none';
axis([0 1 0 1 -0.001 0.004]); title('Numerical Solution'); colorbar;
subplot(3,1,2)
s22 = surf(x,y,s2mat,'FaceAlpha',0.7); s22.EdgeColor = 'none';
axis([0 1 0 1 -0.001 0.004]); title('Exact Solution'); colorbar
subplot(3,1,3)
s23 = surf(x,y,u2mat-s2mat,'FaceAlpha',0.7); s23.EdgeColor = 'none';
title('Numerical - Exact'); colorbar;