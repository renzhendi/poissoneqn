% This function solves Au=f and returns u, run time, and number of iterations.
%
% solverIndex:
% 0 - backslash '\'
% 1 - Jacobi and relaxed Jacobi
% 2 - Gauss-Seidel and SOR
% 3 - SSOR
%
% timingBoolean:
% 0 - no timer
% 1 - launch timer
%
% initGuessType
% 0 - zeros
% 1 - ones
% 2 - rand

function [u,t,iter,errs] = solvers(A, f, solverIndex, timingBoolean, initGuessType, relaxation)

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
uexact = A\f;

if timingBoolean                               % launch timer
    tvec = zeros(1,10);                        % 10 time slots
    switch solverIndex
        case 0                                 % backslash '\'
            tic                                % start timing
            for i = 1:10                       % 10 trials
                u = A\f;
                tvec(i) = toc;                 % cumulative time at the end of each iter
            end
            iter = 1;                          % '\' is not iterative
            errs = 0;                          % '\' yields accurate u
        case 1                                 % Jacobi
            tic;
            for i = 1:10
                [u,iter,errs] = jacobi(A, f, u0, uexact, relaxation);
                tvec(i) = toc;
            end
        case 2
            error('G-S/SOR to be implemented');
        case 3
            error('SSOR to be implemented');
    end
    tvec(2:10) = tvec(2:10) - tvec(1:9);       % break down cumulative time
    t = (sum(tvec) - max(tvec) - min(tvec))/8; % avg time excluding max and min
else                                           % otherwise no timing
    t = nan;
    switch solverIndex
        case 0
            u = A\f;
            iter = 1;
            errs = 0;
        case 1
            [u,iter,errs] = jacobi(A, f, u0, uexact, relaxation);
        case 2
            error('G-S/SOR to be implemented');
        case 3
            error('SSOR to be implemented');
    end
end