% This function constructs and stores the sparse matrix A.

function A = createA(M, N)

h1 = 1/N;                                                % delta x, interval size
h2 = 1/M;                                                % delta y, interval size
diag0 = ones(N-1,1)*2*(h1^2+h2^2)/(h1^2*h2^2);           % trivial step in constructing B
diag1 = -ones(N-1,1)/(h1^2);
B = spdiags([diag1, diag0, diag1], -1:1, N-1, N-1);      % block B, (N-1)*(N-1) matrix
C = -speye(N-1)/(h2^2);                                  % block C, (N-1)*(N-1) matrix
mat_cells = {B,C,zeros(N-1)};                            % trivial step in constructing A
A = cell2mat(mat_cells(toeplitz([1,2,ones(1,M-3)*3])));  % mat A, (M-1)(N-1)*(M-1)(N-1) matrix
save(sprintf('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A/matA_m%d_n%d.mat',M,N),'A');