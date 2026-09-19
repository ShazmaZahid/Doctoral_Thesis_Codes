%% ======================================================================
%  SECOND STAR-METHOD  (GMRES)           Pozza & Zahid, PAMM 24 (2024)
%  "A new Legendre polynomial approach for computing the matrix
%   exponential action on a vector", e202400049.  Section 4.2, Table 2.
%
%  Computes exp(A)*v by solving the full star-linear system
%     (I_MN - (1/2) A (x) T_M) x = v (x) phi_M(-1)
%  iteratively with GMRES, using a matrix-free operator so that the
%  Kronecker matrix is never assembled.
%
%  Requires on the MATLAB path:
%     genCoeffMatrix.m,  and chebfun (used by genCoeffMatrix)
%
%  Uncomment ONE example block below.
%% ======================================================================
clear all;

%% California network
% load California.mat
% A = Problem.A;
% N = length(A);
% b = ones(N,1);
% b = b/norm(b);
% M = 30;
% niter = 100;
% tol = 1e-14;

%% Random Matrix
% N = 500;  % or 1000
% A = rand(N);
% A = A/norm(A);
% N = length(A);
% b = rand(N,1);
% b = b/norm(b);
% M = 10;
% niter = 100;
% tol = 1e-14;

%% Poisson example
J = 50;
[A,b] = poisson_example(J);
b = b/norm(b);
M = 200;
N = length(A);
niter = 300;
tol = 1e-14;

%%
nA = normest(A)

%% Heaviside coefficient matrix
heav = @(t) 1 + 0*t; %defining function
T = genCoeffMatrix(heav,M);  %generating coefficient matrix
T = sparse(T);   %converting coefficient matrix into sparse matrix

% Weight vector (basis)
w = [];
for i=1:M
    w(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2);   %Legendre polynomial basis
end
% truncation
T(end-2:end,:) = sparse(3,M); %to remove error

%% Second star-method: GMRES on the full star-linear system
Atimes = @(v) v - reshape((T/2*reshape(v,M,N))*A.',M*N,1);

tic
[x,flag,relres,iter,resvec] = gmres(Atimes,kron(b,w),[],tol,niter,[],[],kron(b,w));
Star_Sol = sqrt(2)*x(1:M:end);
t_gmres = toc

%% Results (Table 2)
exact=expm(A)*b;

err_star = norm(exact-Star_Sol)/norm(exact)
iter

semilogy(resvec,'bo')
xlabel('GMRES iteration'); ylabel('estimated residual norm');

%% ----------------------------------------------------------------------
function [L,v] = poisson_example(J)
% TREFETHEN, WEIDEMAN, and SCHMELZER 2006
h = 2/J; s = (-1+h:h:1-h)'; % 1D grid
[xx,yy] = meshgrid(s,s); % 2D grid
x = xx(:); y = yy(:); % 2D grid stretched to 1D
L = -gallery('poisson',J-1)/h^2; % 2D Laplacian
I = speye((J-1)^2); % identity
N = 32; theta = pi*(1:2:N-1)/N; % quad pts in (0,pi)
z = N*(.1309-.1194*theta.^2+.2500i*theta); % quad pts on contour
w = N*(-.1194*2*theta+.2500i); % derivatives
c = (1i/N)*exp(z).*w; % quadrature weights
u = (1-x.^2).*(1-y.^2).*exp(x); % initial condition
v = zeros(size(u));
for k = 1:N/2 % quadrature via
    v = v - c(k)*((z(k)*I-0.02*L)\u); % sparse linear solves
end
end
