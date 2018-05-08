function u = vcycle(A, b, u0)

if mod(length(b),2) == 0
    M = diag(diag(A));
    N = M - A;
    us = M\(N*u0 + b);   % pre-smooth
    
    I = speye(size(A,1)/2);
    R = 0.5*kron(I,[1 1]);
    P = 2*R';
    
    Ac = R*A*P;
    rs = b - A*u0;
    rc = R*rs;
    
    ec = vcycle(Ac, rc, zeros(size(rc)));

    es = P*ec;
    u = us+es;
    u = M\(N*u + b);    % post-smooth
else
    u = A\b;
end