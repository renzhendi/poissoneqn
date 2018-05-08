clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script compares 7 iterative methods (Jacobi, relaxed Jacobi, G-S, SOR,
% SSOR, CG, multigrid) when solving the discretized Poisson's equation Au=f.
%
% 1 - Jacobi; 2 - rJacobi(2/3);
% 3 - G-S; 4 - SOR(1.5); 5 - SSOR(1.5);
% 6 - CG; 7 - multigrid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
timingBoolean = 0;
savePlotBoolean = 0;
tol = 10^(-8);

f1=@(x,y)13*pi^2*sin(2*pi*x).*sin(3*pi*y);
f2=@(x,y)-(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5;

solverIndex_list = [1,2,3,4,5,6,7];
N = 66;
M = N;
initGuessType_list = [0];%,1,2];
relaxation_list = [2/3, 1, 1.5]; % Jacobi 2/3 1; SOR/SSOR 1 1.5

%%%%%%%%%%%%%%%%%%%%%%
% RHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
x = (1:N-1)/N;
y = ((1:M-1)/M)';
f1mat = f1(x,y)';           % (N-1)*(M-1) matrix
f1vec = f1mat(:);           % (N-1)(M-1) vector
f2mat = f2(x,y)';
f2vec = f2mat(:);

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
A_name = sprintf('matA_m%d_n%d.mat',M,N);
if exist(A_name,'file') == 2
    load(A_name);
else
    A = createA(M,N);
end
for i = 1:length(initGuessType_list)        % for loop 1: initial condition
    initGuessType = initGuessType_list(i);
    p1 = 10+i;
    p2 = 20+i;
    figure(p1)
    title(sprintf('Problem 1 (initial guess u%d)',initGuessType));
    xlim([0 5000]); %ylim([0 0.1]);
    xlabel('iteration'); ylabel('norm(res)');
    hold on
    figure(p2)
    title(sprintf('Problem 2 (initial guess u%d)',initGuessType));
    xlim([0 5000]); %ylim([0 0.01]);
    xlabel('iteration'); ylabel('norm(res)');
    hold on
    % %s, iter=%d, time=%0.4f (algName,M,N,iter1,t1)
    for k = 1:length(solverIndex_list)      % for loop 2: iterative method
        solverIndex = solverIndex_list(k);
        switch solverIndex
            case 1
                algName = 'Jacobi';
                [u1,t1,iter1,res1] = solvers(A, f1vec, 1, timingBoolean, initGuessType, 1, tol);
                [u2,t2,iter2,res2] = solvers(A, f2vec, 1, timingBoolean, initGuessType, 1, tol);
                figure(p1)
                plot(log(res1),'b');
                figure(p2)
                plot(log(res2),'b');
            case 2
                algName = 'rJacobi(2/3)';
                [u1,t1,iter1,res1] = solvers(A, f1vec, 1, timingBoolean, initGuessType, 2/3, tol);
                [u2,t2,iter2,res2] = solvers(A, f2vec, 1, timingBoolean, initGuessType, 2/3, tol);
                figure(p1)
                plot(log(res1),'r');
                figure(p2)
                plot(log(res2),'r');
            case 3
                tic
                algName = 'GaussSeidel';
                [u1,t1,iter1,res1] = solvers(A, f1vec, 2, timingBoolean, initGuessType, 1, tol);
                [u2,t2,iter2,res2] = solvers(A, f2vec, 2, timingBoolean, initGuessType, 1, tol);
                figure(p1)
                plot(log(res1),'y');
                figure(p2)
                plot(log(res2),'y');
            case 4
                algName = 'SOR(1.5)';
                [u1,t1,iter1,res1] = solvers(A, f1vec, 2, timingBoolean, initGuessType, 1.5, tol);
                [u2,t2,iter2,res2] = solvers(A, f2vec, 2, timingBoolean, initGuessType, 1.5, tol);
                figure(p1)
                plot(log(res1));
                figure(p2)
                plot(log(res2));
            case 5
                algName = 'SSOR(1.5)';
                [u1,t1,iter1,res1] = solvers(A, f1vec, 3, timingBoolean, initGuessType, 1.5, tol);
                [u2,t2,iter2,res2] = solvers(A, f2vec, 3, timingBoolean, initGuessType, 1.5, tol);
                figure(p1)
                plot(log(res1));
                figure(p2)
                plot(log(res2));
            case 6
                algName = 'CG';
                [u1,t1,iter1,res1] = solvers(A, f1vec, 4, timingBoolean, initGuessType, 1, tol);
                [u2,t2,iter2,res2] = solvers(A, f2vec, 4, timingBoolean, initGuessType, 1, tol);
                figure(p1)
                plot(log(res1));
                figure(p2)
                plot(log(res2));
            case 7
                algName = 'multigrid';
                [u1,t1,iter1,res1] = solvers(A, f1vec, 5, timingBoolean, initGuessType, 1, tol);
                [u2,t2,iter2,res2] = solvers(A, f2vec, 5, timingBoolean, initGuessType, 1, tol);
                figure(p1)
                plot(log(res1));
                figure(p2)
                plot(log(res2));
        end
    end
end