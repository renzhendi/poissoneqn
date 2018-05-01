clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script solves general Poisson's equations numerically. It allows users
% to specify M and N (size of the grid), and input f (the RHS of Poisson's
% equation), Then Au=f is solved using the simple finite difference scheme.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 65;                      % number of intervals in x direction
M = 65;                      % number of intervals in y direction
solverIndex = 5;
timingBoolean = 0;
initGuessType = 0;
relaxation = 1;
tol = 10^(-8);
savePlotBoolean = 0;
switch solverIndex
    case 0
        algName = 'backslash';
    case 1
        if relaxation == 1
            algName = 'Jacobi';
        else
            algName = sprintf('rJCB(%0.1f)',relaxation);
        end
    case 2
        if relaxation == 1
            algName = 'GaussSeidel';
        else
            algName = sprintf('SOR(%0.1f)',relaxation);
        end
    case 3
        algName = sprintf('SSOR(%0.1f)',relaxation);
    case 4
        algName = 'CG';
    case 5
        algName = 'Multigrid';
end

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
A_name = sprintf('matA_m%d_n%d.mat',M,N);
if exist(A_name,'file') == 2
    load(A_name);
else
    A = createA(M,N);        % (M-1)(N-1)*(M-1)(N-1) matrix
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
if solverIndex < 4         % splitting methods
    [u1,t1,iter1,res1] = solvers(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);   % (M-1)(N-1) vector and avg run time
    u1mat = vec2mat(u1,N-1);                                                                               % (M-1)*(N-1) matrix
    [u2,t2,iter2,res2] = solvers(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
    u2mat = vec2mat(u2,N-1);
elseif solverIndex == 4    % conjugate gradient
    x0 = zeros(size(f1vec));
    tic
    [u1,iter1,res1] = cg(A, f1vec, x0, A\f1vec, tol);   % (M-1)(N-1) vector and avg run time
    t1 = toc;
    u1mat = vec2mat(u1,N-1);                             % (M-1)*(N-1) matrix
    tic
    [u2,iter2,res2] = cg(A, f2vec, x0, A\f2vec, tol);
    t2 = toc;
    u2mat = vec2mat(u2,N-1);
elseif solverIndex == 5    % multi-grid
    x0 = zeros(size(f1vec));
    tic
    [u1,iter1,res1] = multigrid(A, f1vec, x0, tol);   % (M-1)(N-1) vector and avg run time
    t1 = toc;
    u1mat = vec2mat(u1,N-1);                             % (M-1)*(N-1) matrix
    tic
    [u2,iter2,res2] = multigrid(A, f2vec, x0, tol);
    t2 = toc;
    u2mat = vec2mat(u2,N-1);
else
    error('SolverIndex cannot exceed 5.');
end

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
subplot(3,1,1)        % plot the numerical solution
s11 = surf(x,y,u1mat,'FaceAlpha',0.7); s11.EdgeColor = 'none';
title(sprintf('Eqn1 [%s]: M=%d, N=%d, iter=%d, time=%0.4f',algName,M,N,iter1,t1));
axis([0 1 0 1 -1 1]); colorbar;
xlabel('x'); ylabel('y'); zlabel('numerical soln');
subplot(3,1,2)        % plot the exact solution
s12 = surf(x,y,s1mat,'FaceAlpha',0.7); s12.EdgeColor = 'none';
axis([0 1 0 1 -1 1]); colorbar;
xlabel('x'); ylabel('y'); zlabel('exact soln');
subplot(3,1,3)        % plot the residual (numerical - exact)
s13 = surf(x,y,u1mat-s1mat,'FaceAlpha',0.7); s13.EdgeColor = 'none';
colorbar; xlabel('x'); ylabel('y'); zlabel('numerical - exact');
if savePlotBoolean
    print(sprintf('plots/Eqn1_%s_m%d_n%d.png',algName,M,N),'-dpng');
end

figure(2)
subplot(3,1,1)
s21 = surf(x,y,u2mat,'FaceAlpha',0.7); s21.EdgeColor = 'none';
title(sprintf('Eqn2 [%s]: M=%d, N=%d, iter=%d, time=%0.4f',algName,M,N,iter2,t2));
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
    print(sprintf('plots/Eqn2_%s_m%d_n%d.png',algName,M, N),'-dpng');
end