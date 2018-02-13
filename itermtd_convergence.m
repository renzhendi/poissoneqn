clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script compares 5 iterative methods (Jacobi, relaxed Jacobi, G-S, SOR,
% SSOR) when solving the discretized Poisson's equation system Au=f.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
solverIndex = 1;             % input 1, 2, 3 only
% timingBoolean is 0 by default
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
end

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
N = 50;                      % number of intervals in x direction
M = 50;                      % number of intervals in y direction
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
[~,t1,iter1,errs1,M_mat,N_mat] = solvers(A, f1vec, solverIndex, 0, initGuessType, relaxation, tol);
iter1_vec = 1:iter1;
[~,t2,iter2,errs2] = solvers(A, f2vec, solverIndex, 0, initGuessType, relaxation, tol);
iter2_vec = 1:iter2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical error bound %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
evals = abs(eigs(M_mat\N_mat));
lambda = max(evals);
errs1_bound = log(lambda)*iter1_vec;
errs2_bound = log(lambda)*iter2_vec;

%%%%%%%%%%%%%%%%
% verification %
%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
plot(iter1_vec,log(errs1));
hold on
plot(iter1_vec,errs1_bound)
title(sprintf('Eqn1 [%s, u%d]: M=%d, N=%d, iter=%d, sr=%0.4f',algName,initGuessType,M,N,iter1,lambda));
xlabel('iteration'); ylabel('log(norm(error))');
subplot(2,1,2)
plot(iter2_vec,log(errs2));
hold on
plot(iter2_vec,errs2_bound)
title(sprintf('Eqn2 [%s, u%d]: M=%d, N=%d, iter=%d, sr=%0.4f',algName,initGuessType,M,N,iter2,lambda));
xlabel('iteration'); ylabel('log(norm(error))');
if savePlotBoolean
    print(sprintf('iter_mtd_plots/convergence_%s_u%d_m%d_n%d.png',algName,initGuessType,M,N),'-dpng');
    close;
end