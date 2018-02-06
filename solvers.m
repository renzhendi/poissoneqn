% This function solves Au=f and returns u and run time if required.
%
% solverIndex:
% 0 - backslash '\'
% 1 - Jacobi
% 2 - Gauss-Seidel
% 3 - Relaxed Jacobi
% 4 - SOR
% 5 - SSOR
%
% timingBoolean:
% 0 - no timer
% 1 - launch timer

function [u,t] = solvers(A, f, solverIndex, timingBoolean)

% check number of arguments
% by default, use '\' without timing
if nargin == 2
    solverIndex = 0;
    timingBoolean = 0;
elseif nargin == 3
    timingBoolean = 0;
end

if timingBoolean                               % launch timer
    tvec = zeros(1,10);                        % 10 time slots
    switch solverIndex
        case 0
            tic                                % start timing
            for i = 1:10                       % 10 trials
                u = A\f;
                tvec(i) = toc;                 % cumulative time at the end of each iter
            end
    end
    tvec(2:10) = tvec(2:10) - tvec(1:9);       % break down cumulative time
    t = (sum(tvec) - max(tvec) - min(tvec))/8; % avg time excluding max and min
else                                           % otherwise no timing
    switch solverIndex
        case 0
            u = A\f;
            t = -1;
    end
end