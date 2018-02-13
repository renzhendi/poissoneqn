clear
rng('shuffle')
addpath('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A');
% This script compares 5 iterative methods (Jacobi, relaxed Jacobi, G-S, SOR,
% SSOR) when solving the discretized Poisson's equation system Au=f.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% general initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%
timingBoolean = 0;
savePlotBoolean = 0;
tol = 10^(-8);

f1=@(x,y)13*pi^2*sin(2*pi*x).*sin(3*pi*y);
f2=@(x,y)-(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5;

solverIndex_list = [0, 1, 2, 3];
N_list = [50 100];
M_list = [50 100];
initGuessType_list = [0, 1, 2];
relaxation_list = [0.75, 0.5, 1, 1.5]; % Jacobi 0.75 0.5 1; SOR/SSOR 0.5 1 1.5

for n = 1:length(N_list)        % for loop 1: N
    N = N_list(n);
    x = (1:N-1)/N;
    for m = 1:length(M_list)        % for loop 2: M
        M = M_list(m);
        %%%%%%%%%%%%%%%%%%%%%%
        % RHS initialization %
        %%%%%%%%%%%%%%%%%%%%%%
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
        for i = 1:length(initGuessType_list)        % for loop 3: initial condition
            initGuessType = initGuessType_list(i);
            for k = 1:length(solverIndex_list)        % for loop 4: iterative method
                solverIndex = solverIndex_list(k);
                switch solverIndex
                    case 0
                        algName = 'backslash';
                        [u1a,t1a,iter1a,errs1a] = solvers(A, f1vec, solverIndex, timingBoolean);   % (M-1)(N-1) vector and avg run time
                        % [u1b,t1b,iter1b,errs1b] = solvers2(A, f1vec, solverIndex, timingBoolean);  % component-wise implementation
                        [u2a,t2a,iter2a,errs2a] = solvers(A, f2vec, solverIndex, timingBoolean);
                        % [u2b,t2b,iter2b,errs2b] = solvers2(A, f2vec, solverIndex, timingBoolean);
                        figure
                        subplot(2,1,1)
                        plot(errs1a,'*');
                        title(sprintf('Eqn1 [%s]: M=%d, N=%d, iter=%d, time=%0.4f',algName,M,N,iter1a,t1a));
                        ylim([0 0.1]); xlabel('iteration'); ylabel('norm(error)');
                        subplot(2,1,2)
                        plot(errs2a,'*');
                        title(sprintf('Eqn2 [%s]: M=%d, N=%d, iter=%d, time=%0.4f',algName,M,N,iter2a,t2a));
                        ylim([0 0.01]); xlabel('iteration'); ylabel('norm(error)');
                        if savePlotBoolean
                            print(sprintf('iter_mtd_plots/%s_m%d_n%d.png',algName,M,N),'-dpng');
                            close;
                        end
                    case 1
                        for l = 1:length(relaxation_list)-1        % for loop 5: relaxation
                            relaxation = relaxation_list(l);
                            if relaxation == 1
                                algName = 'Jacobi';
                            else
                                algName = sprintf('rJCB(%0.2f)',relaxation);
                            end
                            [u1a,t1a,iter1a,errs1a] = solvers(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);   % (M-1)(N-1) vector and avg run time
                            % [u1b,t1b,iter1b,errs1b] = solvers2(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);  % component-wise implementation
                            [u2a,t2a,iter2a,errs2a] = solvers(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
                            % [u2b,t2b,iter2b,errs2b] = solvers2(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
                            figure
                            subplot(2,1,1)
                            plot(errs1a);
                            title(sprintf('Eqn1 [%s, u%d]: M=%d, N=%d, iter=%d, time=%0.4f',algName,initGuessType,M,N,iter1a,t1a));
                            ylim([0 0.1]); xlabel('iteration'); ylabel('norm(error)');
                            subplot(2,1,2)
                            plot(errs2a);
                            title(sprintf('Eqn2 [%s, u%d]: M=%d, N=%d, iter=%d, time=%0.4f',algName,initGuessType,M,N,iter2a,t2a));
                            ylim([0 0.01]); xlabel('iteration'); ylabel('norm(error)');
                            if savePlotBoolean
                                print(sprintf('iter_mtd_plots/%s_u%d_m%d_n%d.png',algName,initGuessType,M,N),'-dpng');
                                close;
                            end
                        end
                    case 2
                        for l = 2:length(relaxation_list)
                            relaxation = relaxation_list(l);
                            if relaxation == 1
                                algName = 'GaussSeidel';
                            else
                                algName = sprintf('SOR(%0.2f)',relaxation);
                            end
                            [u1a,t1a,iter1a,errs1a] = solvers(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);   % (M-1)(N-1) vector and avg run time
                            % [u1b,t1b,iter1b,errs1b] = solvers2(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);  % component-wise implementation
                            [u2a,t2a,iter2a,errs2a] = solvers(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
                            % [u2b,t2b,iter2b,errs2b] = solvers2(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
                            figure
                            subplot(2,1,1)
                            plot(errs1a);
                            title(sprintf('Eqn1 [%s, u%d]: M=%d, N=%d, iter=%d, time=%0.4f',algName,initGuessType,M,N,iter1a,t1a));
                            ylim([0 0.1]); xlabel('iteration'); ylabel('norm(error)');
                            subplot(2,1,2)
                            plot(errs2a);
                            title(sprintf('Eqn2 [%s, u%d]: M=%d, N=%d, iter=%d, time=%0.4f',algName,initGuessType,M,N,iter2a,t2a));
                            ylim([0 0.01]); xlabel('iteration'); ylabel('norm(error)');
                            if savePlotBoolean
                                print(sprintf('iter_mtd_plots/%s_u%d_m%d_n%d.png',algName,initGuessType,M,N),'-dpng');
                                close;
                            end
                        end
                    case 3
                        for l = 2:length(relaxation_list)
                            relaxation = relaxation_list(l);
                            algName = sprintf('SSOR(%0.2f)',relaxation);
                            [u1a,t1a,iter1a,errs1a] = solvers(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);   % (M-1)(N-1) vector and avg run time
                            % [u1b,t1b,iter1b,errs1b] = solvers2(A, f1vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);  % component-wise implementation
                            [u2a,t2a,iter2a,errs2a] = solvers(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
                            % [u2b,t2b,iter2b,errs2b] = solvers2(A, f2vec, solverIndex, timingBoolean, initGuessType, relaxation, tol);
                            figure
                            subplot(2,1,1)
                            plot(errs1a);
                            title(sprintf('Eqn1 [%s, u%d]: M=%d, N=%d, iter=%d, time=%0.4f',algName,initGuessType,M,N,iter1a,t1a));
                            ylim([0 0.1]); xlabel('iteration'); ylabel('norm(error)');
                            subplot(2,1,2)
                            plot(errs2a);
                            title(sprintf('Eqn2 [%s, u%d]: M=%d, N=%d, iter=%d, time=%0.4f',algName,initGuessType,M,N,iter2a,t2a));
                            ylim([0 0.01]); xlabel('iteration'); ylabel('norm(error)');
                            if savePlotBoolean
                                print(sprintf('iter_mtd_plots/%s_u%d_m%d_n%d.png',algName,initGuessType,M,N),'-dpng');
                                close;
                            end
                        end
                end
            end
        end
    end
end