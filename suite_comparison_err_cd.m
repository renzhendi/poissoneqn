clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A_cf');
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
n = 65;
timingBoolean = 1;
solverIndex_list = [1,2,3,4,5,6,7];
initGuessType_list = [0,1,2];
formatLegend = '%s (time=%0.4f, iter=%d)';
tol = 10^(-8);

%%%%%%%%%%%%%%%%%%%%%%
% RHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
f=@(x,y)(x.^2-x).*(2*y-3)-2*(y.^2-y);
x = (1:n-1)/n;
y = x';
fmat = f(x,y)';           % (N-1)*(M-1) matrix
fvec = fmat(:);           % (N-1)(M-1) vector

%%%%%%%%%%%%%%%%%%%%%%
% LHS initialization %
%%%%%%%%%%%%%%%%%%%%%%
A_name = sprintf('matA_m%d_n%d.mat',n,n);
if exist(A_name,'file') == 2
    load(A_name);
else
    A = createA(M,n);
end
[~,t1_bs,iter1_bs,err1_bs] = solvers_err(A, fvec, 0, timingBoolean);
leg1_bs = sprintf(formatLegend,'Backslash    ',t1_bs,iter1_bs);

for i = 1:length(initGuessType_list)        % for loop 1: initial condition
    initGuessType = initGuessType_list(i);
    p1 = 10+initGuessType;
    figure(p1)
    title(sprintf('Convection-diffusion Equation (initial guess u%d)',initGuessType));
    xlim([0 5000]); ylim([-18 6]);
    xlabel('iteration'); ylabel('log[norm(err)]');
    hold on
    for k = 1:length(solverIndex_list)      % for loop 2: iterative method
        solverIndex = solverIndex_list(k);
        switch solverIndex
            case 1 % 'Jacobi';
                [~,t1_jcb,iter1_jcb,err1_jcb] = solvers_err(A, fvec, 1, timingBoolean, initGuessType, 1, tol);
                leg1_jcb = sprintf(formatLegend,'Jacobi          ',t1_jcb,iter1_jcb);
            case 2 % 'rJacobi(2/3)';
                [~,t1_rj,iter1_rj,err1_rj] = solvers_err(A, fvec, 1, timingBoolean, initGuessType, 2/3, tol);
                leg1_rj = sprintf(formatLegend,'rJacobi(2/3) ',t1_rj,iter1_rj);
            case 3 % 'GaussSeidel';
                [~,t1_gs,iter1_gs,err1_gs] = solvers_err(A, fvec, 2, timingBoolean, initGuessType, 1, tol);
                leg1_gs = sprintf(formatLegend,'GaussSeidel',t1_gs,iter1_gs);
            case 4 % 'SOR(3/2)';
                [~,t1_sor,iter1_sor,err1_sor] = solvers_err(A, fvec, 2, timingBoolean, initGuessType, 1.5, tol);
                leg1_sor = sprintf(formatLegend,'SOR(3/2)     ',t1_sor,iter1_sor);
            case 5 % 'SSOR(3/2)';
                [~,t1_ss,iter1_ss,err1_ss] = solvers_err(A, fvec, 3, timingBoolean, initGuessType, 1.5, tol);
                leg1_ss = sprintf(formatLegend,'SSOR(3/2)   ',t1_ss,iter1_ss);
            case 6 % 'CG';
                [~,t1_cg,iter1_cg,err1_cg] = solvers_err(A, fvec, 4, timingBoolean, initGuessType, 1, tol);
                leg1_cg = sprintf(formatLegend,'CG               ',t1_cg,iter1_cg);
            case 7 % 'multigrid';
                [~,t1_mg,iter1_mg,err1_mg] = solvers_err(A, fvec, 5, timingBoolean, initGuessType, 1, tol);
                leg1_mg = sprintf(formatLegend,'Multigrid       ',t1_mg,iter1_mg);
        end
    end
    figure(p1)
    plot(log(err1_bs),'o','LineWidth',2); plot(log(err1_jcb),'LineWidth',2);
    plot(log(err1_rj),'LineWidth',2); plot(log(err1_gs),'LineWidth',2);
    plot(log(err1_sor),'LineWidth',2); plot(log(err1_ss),'LineWidth',2);
    plot(log(err1_cg),'LineWidth',2);plot(log(err1_mg),'LineWidth',2);
    legend(leg1_bs,leg1_jcb,leg1_rj,leg1_gs,leg1_sor,leg1_ss,leg1_cg,leg1_mg);
end