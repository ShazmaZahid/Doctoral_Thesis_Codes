%% ======================================================================
%  FIRST STAR-METHOD                     Pozza & Zahid, PAMM 24 (2024)
%  "A new Legendre polynomial approach for computing the matrix
%   exponential action on a vector", e202400049.  Section 4.1, Table 1.
%
%  Computes exp(A)*v by projecting A onto a Krylov subspace (Arnoldi) and
%  solving the reduced Stein equation with dlyap.  Compared against the
%  classical Krylov approximation Q_n*expm(H_n)*e_1.
%
%  Requires on the MATLAB path:
%     Arnoldimgs.m, genCoeffMatrix.m,  and chebfun (used by genCoeffMatrix)
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
% k = 22;

%% Random Matrix
% N = 500;  % or 1000
% A = rand(N);
% A = A/norm(A);
% N = length(A);
% b = rand(N,1);
% b = b/norm(b);
% M = 10;
% k = 10;   % 6, 8, 10 for N=500 ; 4, 6, 8 for N=1000

%% Poisson example
J = 50;
[A,b] = poisson_example(J);
b = b/norm(b);
k = 150;  % Number of Arnoldi iterations
M = 200;
N = length(A);

%%
nA = normest(A)

%% Arnoldi
tic
[Vk,Hk]=Arnoldimgs(A,b,k);
Hk = Hk(1:k-1,1:k-1);
toc

%% Pade Krylov Approximation  (classical baseline)
e1 = zeros(k-1,1);
e1(1) = 1;
tic
Pade_Sol=Vk(:,1:k-1)*expm(Hk)*e1;
toc

%% Star Krylov Approximation  (first star-method)
% Heaviside function discretization
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
[U,S] = schur(full(T/2),'complex');
Uw = (U'*w);

tic
Y = dlyap(S,Hk(1:k-1,1:k-1).', Uw*e1');
X = U*Y;
Star_Sol = sqrt(2)*Vk(:,1:k-1)*X(1,:).';
toc

%% Results (Table 1)
exact=expm(A)*b;

err_pade = norm(exact-Pade_Sol)/norm(exact)
err_star = norm(exact-Star_Sol)/norm(exact)

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
