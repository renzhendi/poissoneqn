% This function constructs and stores the sparse matrix A.

function A = createA_cf(n,a1,a2)

h = 1/n;                                                 % interval size
diag0 = ones(N-1,1)*4/(h^2);                             % trivial step in constructing B
diag_l = -ones(N-1,1)*(1+a1*h/2)/(h^2);
diag_u = -ones(N-1,1)*(1-a1*h/2)/(h^2);
B = spdiags([diag_l, diag0, diag_u], -1:1, N-1, N-1);    % block B, (N-1)*(N-1) matrix
C_l = -speye(N-1)*(1+a2*h/2)/(h2^2);                     % block C, (N-1)*(N-1) matrix
C_u = -speye(N-1)*(1-a2*h/2)/(h2^2);
mat_cells = {B,C_l,C_u,zeros(N-1)};                                          % trivial step in constructing A
A = cell2mat(mat_cells(toeplitz([1,2,ones(1,n-3)*4],[1,3,ones(1,n-3)*4])));  % mat A, (M-1)(N-1)*(M-1)(N-1) matrix
save(sprintf('C:/Users/Œ‚ﬁ»ïF/Documents/MATLAB/MMSC/poissoneqn/A_cf/matA_m%d_n%d.mat',M,N),'A');