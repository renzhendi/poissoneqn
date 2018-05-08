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
N = 65;
M = 65;
timingBoolean = 1;
solverIndex_list = [1,2,3,4,5,6,7];
initGuessType_list = [0,1,2];
formatLegend = '%s (time=%0.4f, iter=%d)';
tol = 10^(-8);

%%%%%%%%%%%%%%%%%%%%%%
% RHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
f1=@(x,y)13*pi^2*sin(2*pi*x).*sin(3*pi*y);
f2=@(x,y)-(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5;
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
[~,t1_bs,iter1_bs,err1_bs] = solvers_err(A, f1vec, 0, timingBoolean);
[~,t2_bs,iter2_bs,err2_bs] = solvers_err(A, f2vec, 0, timingBoolean);
leg1_bs = sprintf(formatLegend,'Backslash    ',t1_bs,iter1_bs);
leg2_bs = sprintf(formatLegend,'Backslash    ',t2_bs,iter2_bs);

for i = 1:length(initGuessType_list)        % for loop 1: initial condition
    initGuessType = initGuessType_list(i);
    p1 = 10+initGuessType;
    p2 = 20+initGuessType;
    figure(p1)
    title(sprintf('Problem 1 (initial guess u%d)',initGuessType));
    xlim([0 5000]); ylim([-18 6]);
    xlabel('iteration'); ylabel('log[norm(err)]');
    hold on
    figure(p2)
    title(sprintf('Problem 2 (initial guess u%d)',initGuessType));
    xlim([0 5000]); ylim([-18 6]);
    xlabel('iteration'); ylabel('log[norm(err)]');
    hold on
    for k = 1:length(solverIndex_list)      % for loop 2: iterative method
        solverIndex = solverIndex_list(k);
        switch solverIndex
            case 1 % 'Jacobi';
                [~,t1_jcb,iter1_jcb,err1_jcb] = solvers_err(A, f1vec, 1, timingBoolean, initGuessType, 1, tol);
                [~,t2_jcb,iter2_jcb,err2_jcb] = solvers_err(A, f2vec, 1, timingBoolean, initGuessType, 1, tol);
                leg1_jcb = sprintf(formatLegend,'Jacobi          ',t1_jcb,iter1_jcb);
                leg2_jcb = sprintf(formatLegend,'Jacobi          ',t2_jcb,iter2_jcb);
            case 2 % 'rJacobi(2/3)';
                [~,t1_rj,iter1_rj,err1_rj] = solvers_err(A, f1vec, 1, timingBoolean, initGuessType, 2/3, tol);
                [~,t2_rj,iter2_rj,err2_rj] = solvers_err(A, f2vec, 1, timingBoolean, initGuessType, 2/3, tol);
                leg1_rj = sprintf(formatLegend,'rJacobi(2/3) ',t1_rj,iter1_rj);
                leg2_rj = sprintf(formatLegend,'rJacobi(2/3) ',t2_rj,iter2_rj);
            case 3 % 'GaussSeidel';
                [~,t1_gs,iter1_gs,err1_gs] = solvers_err(A, f1vec, 2, timingBoolean, initGuessType, 1, tol);
                [~,t2_gs,iter2_gs,err2_gs] = solvers_err(A, f2vec, 2, timingBoolean, initGuessType, 1, tol);
                leg1_gs = sprintf(formatLegend,'GaussSeidel',t1_gs,iter1_gs);
                leg2_gs = sprintf(formatLegend,'GaussSeidel',t2_gs,iter2_gs);
            case 4 % 'SOR(3/2)';
                [~,t1_sor,iter1_sor,err1_sor] = solvers_err(A, f1vec, 2, timingBoolean, initGuessType, 1.5, tol);
                [~,t2_sor,iter2_sor,err2_sor] = solvers_err(A, f2vec, 2, timingBoolean, initGuessType, 1.5, tol);
                leg1_sor = sprintf(formatLegend,'SOR(3/2)     ',t1_sor,iter1_sor);
                leg2_sor = sprintf(formatLegend,'SOR(3/2)     ',t2_sor,iter2_sor);
            case 5 % 'SSOR(3/2)';
                [~,t1_ss,iter1_ss,err1_ss] = solvers_err(A, f1vec, 3, timingBoolean, initGuessType, 1.5, tol);
                [~,t2_ss,iter2_ss,err2_ss] = solvers_err(A, f2vec, 3, timingBoolean, initGuessType, 1.5, tol);
                leg1_ss = sprintf(formatLegend,'SSOR(3/2)   ',t1_ss,iter1_ss);
                leg2_ss = sprintf(formatLegend,'SSOR(3/2)   ',t2_ss,iter2_ss);
            case 6 % 'CG';
                [~,t1_cg,iter1_cg,err1_cg] = solvers_err(A, f1vec, 4, timingBoolean, initGuessType, 1, tol);
                [~,t2_cg,iter2_cg,err2_cg] = solvers_err(A, f2vec, 4, timingBoolean, initGuessType, 1, tol);
                leg1_cg = sprintf(formatLegend,'CG               ',t1_cg,iter1_cg);
                leg2_cg = sprintf(formatLegend,'CG               ',t2_cg,iter2_cg);
            case 7 % 'multigrid';
                [~,t1_mg,iter1_mg,err1_mg] = solvers_err(A, f1vec, 5, timingBoolean, initGuessType, 1, tol);
                [~,t2_mg,iter2_mg,err2_mg] = solvers_err(A, f2vec, 5, timingBoolean, initGuessType, 1, tol);
                leg1_mg = sprintf(formatLegend,'Multigrid       ',t1_mg,iter1_mg);
                leg2_mg = sprintf(formatLegend,'Multigrid       ',t2_mg,iter2_mg);
        end
    end
    figure(p1)
    plot(log(err1_bs),'o','LineWidth',2); plot(log(err1_jcb),'LineWidth',2);
    plot(log(err1_rj),'LineWidth',2); plot(log(err1_gs),'LineWidth',2);
    plot(log(err1_sor),'LineWidth',2); plot(log(err1_ss),'LineWidth',2);
    plot(log(err1_cg),'LineWidth',2); plot(log(err1_mg),'LineWidth',2);
    legend(leg1_bs,leg1_jcb,leg1_rj,leg1_gs,leg1_sor,leg1_ss,leg1_cg,leg1_mg);
    
    figure(p2)
    plot(log(err2_bs),'o','LineWidth',2); plot(log(err2_jcb),'LineWidth',2);
    plot(log(err2_rj),'LineWidth',2); plot(log(err2_gs),'LineWidth',2);
    plot(log(err2_sor),'LineWidth',2); plot(log(err2_ss),'LineWidth',2);
    plot(log(err2_cg),'LineWidth',2); plot(log(err2_mg),'LineWidth',2);
    legend(leg2_bs,leg2_jcb,leg2_rj,leg2_gs,leg2_sor,leg2_ss,leg2_cg,leg2_mg);
end