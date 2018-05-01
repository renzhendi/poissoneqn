clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script compares 5 iterative methods (Jacobi, relaxed Jacobi, G-S, SOR,
% SSOR) when solving the discretized Poisson's equation system Au=f.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
solverIndex = 4;             % input 1, 2, 3 only
% timingBoolean is 0 by default
initGuessType = 0;
relaxation = 1.5;
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
N = 200;                      % number of intervals in x direction
M = 200;                      % number of intervals in y direction
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
x0 = zeros(size(f1vec));
[u1,iter1,errs1,errsA1,res1] = cg_reserr(A, f1vec, x0, A\f1vec, tol);   % (M-1)(N-1) vector and avg run time
iter1_vec = 0:iter1-1;                                                  % recorded u0 in iter1
[u2,iter2,errs2,errsA2,res2] = cg_reserr(A, f2vec, x0, A\f2vec, tol);
iter2_vec = 0:iter2-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical error bound %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
evals = abs(eig(A));
lambda1 = max(evals);
lambda2 = min(evals);
k = lambda1/lambda2;
errs1_bound = 2*errsA1(1)*((sqrt(k)-1)/(sqrt(k)+1)).^iter1_vec;
errs2_bound = 2*errsA2(1)*((sqrt(k)-1)/(sqrt(k)+1)).^iter2_vec;

%%%%%%%%%%%%%%%%
% verification %
%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
plot(iter1_vec,log(errsA1(1:end)));
hold on
plot(iter1_vec,log(errs1_bound))
title(sprintf('Eqn1 [%s, u%d]: M=%d, N=%d, iter=%d, kappa=%0.4f',algName,initGuessType,M,N,iter1,k));
xlabel('iteration'); ylabel('A-norm(error)');
subplot(2,1,2)
plot(iter2_vec,log(errsA2(1:end)));
hold on
plot(iter2_vec,log(errs2_bound));
title(sprintf('Eqn2 [%s, u%d]: M=%d, N=%d, iter=%d, kappa=%0.4f',algName,initGuessType,M,N,iter2,k));
xlabel('iteration'); ylabel('A-norm(error)');
if savePlotBoolean
    print(sprintf('iter_mtd_plots/convergence_%s_u%d_m%d_n%d.png',algName,initGuessType,M,N),'-dpng');
    close;
end