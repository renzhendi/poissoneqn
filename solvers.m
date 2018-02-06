% This function solves Au=f and returns u and run time if required.

function [u,t] = solvers(A, f, solverIndex, timingBoolean)

% check number of arguments
% by default, using '\' and no timing
if nargin == 2
    solverIndex = 1;
    timingBoolean = 0;
elseif nargin == 3
    timingBoolean = 0;
end

if timingBoolean
    switch solverIndex
        case 1
            tic
            for i = 1:10
                u = A\f;
            end
            t = toc/10;
    end
    % otherwise no timing
else
    switch solverIndex
        case 1
            u = A\f;
            t = -1;
    end
end