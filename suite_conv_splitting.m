clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script compares 5 iterative methods (Jacobi, relaxed Jacobi, G-S, SOR,
% SSOR) when solving the discretized Poisson's equation system Au=f.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 65;                      % number of intervals in x direction
M = 65;                      % number of intervals in y direction
solverIndex = 1;             % input 1, 2, 3 only
initGuessType = 2;
relaxation = 1;
tol = 10^(-8);

switch solverIndex
    case 0
        algName = 'backslash';
    case 1
        if relaxation == 1
            algName = 'Jacobi';
        else
            algName = 'rJacobi(2/3)';
        end
    case 2
        if relaxation == 1
            algName = 'GaussSeidel';
        else
            algName = 'SOR(3/2)';
        end
    case 3
        algName = 'SSOR(3/2)';
end

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
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
[~,t1,iter1,errs1,M_mat,N_mat] = solvers_err(A, f1vec, solverIndex, 0, initGuessType, relaxation, tol);
iter1_vec = 1:iter1;
[~,t2,iter2,errs2] = solvers_err(A, f2vec, solverIndex, 0, initGuessType, relaxation, tol);
iter2_vec = 1:iter2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical error bound %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
evals = abs(eigs(M_mat\N_mat));
lambda = max(evals);
errs1_bound = log(lambda)*iter1_vec + max(0,log(errs1(1)));
errs2_bound = log(lambda)*iter2_vec + max(0,log(errs2(1)));

%%%%%%%%%%%%%%%%
% verification %
%%%%%%%%%%%%%%%%
figure(1)
plot(iter1_vec,log(errs1),'LineWidth',2);
hold on
plot(iter1_vec,errs1_bound,'LineWidth',2);
title(sprintf('Problem 1 (%s, u%d, lambda=%0.4f)',algName,initGuessType,lambda));
xlabel('iteration'); ylabel('log[norm(err)]');
legend('Measured error','Theoretical bound');

figure(2)
plot(iter2_vec,log(errs2),'LineWidth',2);
hold on
plot(iter2_vec,errs2_bound,'LineWidth',2);
title(sprintf('Problem 2 (%s, u%d, lambda=%0.4f)',algName,initGuessType,lambda));
xlabel('iteration'); ylabel('log[norm(err)]');
legend('Measured error','Theoretical bound');