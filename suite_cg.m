clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script compares 5 iterative methods (Jacobi, relaxed Jacobi, G-S, SOR,
% SSOR) when solving the discretized Poisson's equation system Au=f.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
algName = 'CG';
initGuessType = 0;
tol = 10^(-8);
savePlotBoolean = 0;

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
N = 100;                      % number of intervals in x direction
M = 100;                      % number of intervals in y direction
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
[u1,iter1,errA1,res1] = cg_A(A,f1vec,x0,tol);   % (M-1)(N-1) vector and avg run time
iter1_vec = 0:iter1-1;                          % recorded u0 in iter1
[u2,iter2,errA2,res2] = cg_A(A,f2vec,x0,tol);
iter2_vec = 0:iter2-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical error bound %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
evals = abs(eig(A));
lambda1 = max(evals);
lambda2 = min(evals);
k = lambda1/lambda2;
errs1_bound = 2*errA1(1)*((sqrt(k)-1)/(sqrt(k)+1)).^iter1_vec;
errs2_bound = 2*errA2(1)*((sqrt(k)-1)/(sqrt(k)+1)).^iter2_vec;

%%%%%%%%%%%%%%%%
% verification %
%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
plot(iter1_vec,log(errA1(1:end)));
hold on
plot(iter1_vec,log(errs1_bound))
title(sprintf('Eqn1 [%s, u%d]: M=%d, N=%d, iter=%d, kappa=%0.4f',algName,initGuessType,M,N,iter1,k));
xlabel('iteration'); ylabel('A-norm(error)');
subplot(2,1,2)
plot(iter2_vec,log(errA2(1:end)));
hold on
plot(iter2_vec,log(errs2_bound));
title(sprintf('Eqn2 [%s, u%d]: M=%d, N=%d, iter=%d, kappa=%0.4f',algName,initGuessType,M,N,iter2,k));
xlabel('iteration'); ylabel('A-norm(error)');
if savePlotBoolean
    print(sprintf('iter_mtd_plots/convergence_%s_u%d_m%d_n%d.png',algName,initGuessType,M,N),'-dpng');
    close;
end