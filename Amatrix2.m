function A = Amatrix2(N,M)
hx=1/N;
hy=1/M;

C=-1/hy^2*eye(N-1,N-1);
B=toeplitz([2/hx^2+2/hy^2,-1/hx^2,zeros(1,N-3)]);
A=kron(eye(M-1),B)+kron(toeplitz([0,1,zeros(1,M-3)]),C);

% A=zeros((n-1)*(m-1),(n-1)*(m-1));
% C=-1/hy^2*eye(n-1,n-1);
% B=toeplitz([2/hx^2+2/hy^2,-1/hx^2,zeros(1,n-3)]);
% 
% for i=1:m-1
%     A((i-1)*(n-1)+1:(n-1)*i,(i-1)*(n-1)+1:(n-1)*i)=B;
%     if (i ~= m-1) A((i-1)*(n-1)+1:(n-1)*i,(i)*(n-1)+1:(n-1)*(i+1))=C; end %upper diag 
%     if (i ~= m-1) A((i)*(n-1)+1:(n-1)*(i+1),(i-1)*(n-1)+1:(n-1)*i)=C; end %lower diag
% end