clear
% u_{xx} + u_{yy} = f(x,y)
% on [a,b] x [a,b].  


a = 0; 
b = 1; 
N = 20; 
M = 20;
h_x = (b-a)/N;
h_y = (b-a)/M;
x = linspace(a,b,N+1);   % grid points x including boundaries
y = linspace(a,b,M+1);   % grid points y including boundaries


[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

i = 2:N;            % indices of interior points in x
j = 2:M;            % indices of interior points in y
Xin = X(i,j);       % interior points
Yin = Y(i,j);

%f=@(x,y) 13*pi^2.*sin(2*pi*x).*sin(3*pi*y);   %forced function for right hand side
f=@(x,y) -1*(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5;
%-(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5

%u_exact = @(x,y) sin(2*pi*x).*sin(3*pi*y);         % exact solution of u
u_exact = @(x,y) (x-1).^5.*x.^2.*y.*(y-1); 

rhs = f(Xin,Yin);        % evaluate f at interior points for right hand side
                           % rhs is modified below for boundary conditions.
%real solution
ureal= u_exact(Xin,Yin); 
Ureal=zeros(N+1,M+1);
Ureal(2:N,2:M)=ureal;
urealcol=reshape(ureal,(N-1)*(M-1),1);

% convert the 2d grid function rhs into a column vector for rhs of system:
F=zeros(N+1,M+1);
F(2:N,2:M)=rhs;
Fcol=reshape(rhs,(N-1)*(M-1),1);

% form matrix A:

ex = ones(N-1,1);
ey = ones(M-1,1);
B = spdiags([-ex/(h_x^2) (2/(h_x^2)+2/(h_y^2))*ex -ex/(h_x^2)],[-1 0 1],N-1,N-1); %inner block B is acting on x
C = spdiags([ey ey],[-1 1],M-1,M-1);
I = speye(M-1,M-1);
I_1 = -speye(N-1,N-1)/(h_y^2);
A = (kron(I,B) + kron(C,I_1));

%extract strict lower triangle; dia; uppertriangle

L=tril(A,-1);
D=diag(diag(A));
R=triu(A,1);

%Gauss Seidel and SOR
w=1;
% initial guess of the interior points
u0=zeros((N-1)*(M-1));
%% 
Uiteration_GS=u0;
k_GS=0;  %kth column corresponding to k-1 iteratin
err_GS = max(abs(Uiteration_GS(:,k_GS+1)-urealcol));

while k_GS<50000 && err_GS > 10e-4
    %(D+WL)u(k+1)=-Ru(k)+f
    unext=(D+w*L)\(((1-w)*D-w*R)*Uiteration_GS(:,k_GS+1)+w*Fcol);
    Uiteration_GS=[Uiteration_GS unext];
    err_GS = norm(Uiteration_GS(:,k_GS+2)-urealcol);
    k_GS=k_GS+1;
end 
    step=k_GS
    u_GSsolution_col=Uiteration_GS(:,k_GS+1);
    uGSgrid=zeros(N+1,M+1);
    uGSgrid(2:N,2:M) = reshape(u_GSsolution_col,N-1,M-1);
    fprintf('Error of GS to true solution of PDE = %10.3e \n',err_GS)

% %     

% %% 
% w=1;
% % initial guess of the interior points
% u0=zeros((N-1)*(M-1));
% Uiteration_SSOR=u0;
% k_SSOR=0;  %kth column corresponding to k-1 iteratin
% err_SSOR = max(abs(Uiteration_SSOR(:,k_SSOR+1)-urealcol));
% 
% while k_SSOR<50000 && err_SSOR > 10e-4
%     %(D+WL)u(k+1)=-Ru(k)+f
%     u_half=(D+w*L)\(-((1-w)*D-w*R)*Uiteration_SSOR(:,k_SSOR+1)+w*Fcol);
%     u_full=(D+w*R)\(((1-w)*D-w*L)*u_half+w*Fcol);
%     Uiteration_SSOR=[Uiteration_SSOR u_full];
%     err_SSOR = max(abs(Uiteration_SSOR(:,k_SSOR+2)-urealcol));
%     k_SSOR=k_SSOR+1;
% end 
%     step=k_SSOR
%     u_SSORsolution_col=Uiteration_SSOR(:,k_SSOR+1);
%     uSSORgrid=zeros(N+1,M+1);
%     uSSORgrid(2:N,2:M) = reshape(u_SSORsolution_col,N-1,M-1);
%     fprintf('Error of GS to true solution of PDE = %10.3e \n',err_SSOR)
% 
% %err_k= max(abs(uGS_k-urealcol))>10e-4;
% %if err_
%% 



% Solve the linear system:
usolve = A\Fcol;  

% reshape vector solution uvec as a grid function and 
% insert this interior solution into usoln for plotting purposes:
% (recall boundary conditions in usoln are already set) 

usolvegrid=zeros(N+1,M+1);
usolvegrid(2:N,2:M) = reshape(usolve,N-1,M-1);

% assuming true solution is known and stored in utrue:
err = max(abs(usolve-urealcol));
fprintf('Error relative to true solution of PDE = %10.3e \n',err)

% plot results:
%% 

% plot exact solution 
figure()
surf(Ureal)
axis tight
%% 

% plot approx. solution backslash
figure()
surf(usolvegrid)
axis tight
%% 
% plot approx. solution GS
figure()
surf(uGSgrid)
axis tight

% %% 
% % plot approx. solution SSOR
% figure()
% surf(uSSORgrid)
% axis tight
