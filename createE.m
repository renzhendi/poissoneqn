% This function constructs and stores the sparse matrix A.

function evals = createE(A,M,N)

evals = abs(eig(A));
save(sprintf('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A/evalA_m%d_n%d.mat',M,N),'evals');