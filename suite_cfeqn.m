clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A_cf');
% This script solves general Poisson's equations numerically. It allows users
% to specify M and N (size of the grid), and input f (the RHS of Poisson's
% equation), Then Au=f is solved using the simple finite difference scheme.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 65;                      % number of intervals in x&y directions
solverIndex = 0;
timingBoolean = 0;
initGuessType = 0;
relaxation = 1;
tol = 10^(-8);
savePlotBoolean = 0;
switch solverIndex
    case 0
        algName = 'Backslash';
    case 1
        if relaxation == 1
            algName = 'Jacobi';
        else
            algName = 'rJacobi(2/3)';
            % algName = sprintf('rJCB(%0.1f)',relaxation);
        end
    case 2
        if relaxation == 1
            algName = 'GaussSeidel';
        else
            algName = 'SOR(3/2)';
            % algName = sprintf('SOR(%0.1f)',relaxation);
        end
    case 3
        algName = 'SSOR(3/2)';
        % algName = sprintf('SSOR(%0.1f)',relaxation);
    case 4
        algName = 'CG';
    case 5
        algName = 'Multigrid';
end

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
A_name = sprintf('matA_m%d_n%d.mat',n,n);
if exist(A_name,'file') == 2
    load(A_name);
else
    A = createA(n,n);        % (M-1)(N-1)*(M-1)(N-1) matrix
end

%%%%%%%%%%%%%%%%%%%%%%
% RHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
x = (1:n-1)/n;                        % horizontal axis
y = x';                               % vertical axis
f=@(x,y)(x.^2-x).*(2*y-3)-2*(y.^2-y);
fmat = f(x,y)';                       % (n-1)*(n-1) matrix
fvec = fmat(:);                       % (n-1)(n-1) vector

%%%%%%%%%%%%%%%%%%
% numerical soln %
%%%%%%%%%%%%%%%%%%
if solverIndex <= 5
    [u,t,iter,res] = solvers(A, fvec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
    u1mat = vec2mat(u,n-1);
else
    error('SolverIndex cannot exceed 5.');
end

%%%%%%%%%%%%%%
% exact soln %
%%%%%%%%%%%%%%
s=@(x,y)(x.^2-x).*(y.^2-y);
smat = s(x,y);

%%%%%%%%%%%%%%%%
% verification %
%%%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1)        % plot the numerical solution
s1 = surf(x,y,u1mat,'FaceAlpha',0.7); s1.EdgeColor = 'none';
title(sprintf('Convection-diffusion Equation (using %s)',algName));
% title(sprintf('Problem 1 [%s]: M=%d, N=%d, iter=%d, time=%0.4f',algName,M,N,iter1,t1));
axis([0 1 0 1 0 0.065]); colorbar;
xlabel('x'); ylabel('y'); zlabel('numerical soln');
subplot(3,1,2)        % plot the exact solution
s2 = surf(x,y,smat,'FaceAlpha',0.7); s2.EdgeColor = 'none';
axis([0 1 0 1 0 0.065]); colorbar;
xlabel('x'); ylabel('y'); zlabel('exact soln');
subplot(3,1,3)        % plot the residual (numerical - exact)
s3 = surf(x,y,u1mat-smat,'FaceAlpha',0.7); s3.EdgeColor = 'none';
colorbar; xlabel('x'); ylabel('y'); zlabel('numerical - exact');
if savePlotBoolean
    print(sprintf('P1_%s.png',algName),'-dpng');
    % print(sprintf('P1_%s_m%d_n%d.png',algName,M,N),'-dpng');
end