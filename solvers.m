% This function solves Au=f and returns
% u, solution;
% t, run time;
% iter, number of iterations;
% res_vec, vector of residuals;
% M & N, iteration matrices if choosing splitting methods.
%
% solverIndex:
% 0 - backslash '\'
% 1 - Jacobi and relaxed Jacobi
% 2 - Gauss-Seidel and SOR
% 3 - SSOR
% 4 - CG
% 5 - Multigrid
%
% timingBoolean:
% 0 - no timer
% 1 - launch timer
%
% initGuessType
% 0 - zeros
% 1 - ones
% 2 - rand

function [u,t,iter,res_vec,M,N] = solvers(A,f,solverIndex,timingBoolean,initGuessType,relaxation,tol)

if nargin < 7
    tol = 10^(-8);
end
if nargin < 6
    relaxation = 1;     % by default, no relaxation
end
if nargin < 5
    initGuessType = 0;  % by default, u0=0
end
if nargin < 4
    timingBoolean = 0;  % by default, no timing
end
if nargin < 3
    solverIndex = 0;    % by default, use '\'
end
if solverIndex > 0
    switch initGuessType
        case 0          % u0=0
            u0 = 0*f;
        case 1          % u0=1
            u0 = 0*f+1;
        case 2          % u0=rand
            u0 = rand(size(f));
    end
end

if timingBoolean                               % launch timer
    M = -1;                                    % does record iter mat in timing mode
    N = -1;
    tvec = zeros(1,10);                        % 10 time slots
    switch solverIndex
        case 0                                 % backslash '\'
            tic                                % start timing
            for i = 1:10                       % 10 trials
                u = A\f;
                tvec(i) = toc;                 % cumulative time at the end of each iter
            end
            iter = 1;                          % '\' is not iterative
            res_vec = 0;                          % '\' yields accurate u
        case 1                                 % Jacobi
            tic;
            for i = 1:10
                [u,iter,res_vec] = jacobi(A, f, u0, relaxation, tol);
                tvec(i) = toc;
            end
        case 2
            tic;
            for i = 1:10
                [u,iter,res_vec] = gaussseidel(A, f, u0, relaxation, tol);
                tvec(i) = toc;
            end
        case 3
            tic;
            for i = 1:10
                [u,iter,res_vec] = ssor(A, f, u0, relaxation, tol);
                tvec(i) = toc;
            end
        case 4
            tic;
            for i = 1:10
                [u,iter,res_vec] = cg(A, f, u0, tol);
                tvec(i) = toc;
            end
        case 5
            tic;
            for i = 1:10
                [u,iter,res_vec] = multigrid(A, f, u0, tol);
                tvec(i) = toc;
            end
    end
    tvec(2:10) = tvec(2:10) - tvec(1:9);       % break down cumulative time
    t = (sum(tvec) - max(tvec) - min(tvec))/8; % avg time excluding max and min
else                                           % otherwise no timing
    t = nan;
    switch solverIndex
        case 0
            u = A\f;
            iter = 1;
            res_vec = 0;
        case 1
            [u,iter,res_vec,M,N] = jacobi(A, f, u0, relaxation, tol);
        case 2
            [u,iter,res_vec,M,N] = gaussseidel(A, f, u0, relaxation, tol);
        case 3
            [u,iter,res_vec,M,N] = ssor(A, f, u0, relaxation, tol);
        case 4
            [u,iter,res_vec] = cg(A, f, u0, tol);
            M = -1;
            N = -1;
        case 5
            [u,iter,res_vec] = multigrid(A, f, u0, tol);
            M = -1;
            N = -1;
    end
end