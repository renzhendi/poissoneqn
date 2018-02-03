function A = Amatrix(N)
h=1/N;

C=-1/h^2*eye(N-1,N-1);
B=toeplitz([4,-1,zeros(1,N-3)])/1/h^2;
A=kron(eye(N-1),B)+kron(toeplitz([0,1,zeros(1,N-3)]),C);

% A=zeros((n-1)^2,(n-1)^2);
% C=-1/h^2*eye(n-1,n-1);
% B=toeplitz([4,-1,zeros(1,n-3)])/1/h^2;
% 
% for i=1:n-1
%     A((i-1)*(n-1)+1:(n-1)*i,(i-1)*(n-1)+1:(n-1)*i)=B;
%     if (i ~= n-1) A((i-1)*(n-1)+1:(n-1)*i,(i)*(n-1)+1:(n-1)*(i+1))=C; end %upper diag 
%     if (i ~= n-1) A((i)*(n-1)+1:(n-1)*(i+1),(i-1)*(n-1)+1:(n-1)*i)=C; end %lower diag
% end

    
