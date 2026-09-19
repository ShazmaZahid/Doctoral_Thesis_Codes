%% ======================================================================
%  TABLE 2  (PAMM 2024)  --  SECOND STAR-METHOD (GMRES)
%  Runs every case in one go and prints the table.
%
%  Requires on the path: genCoeffMatrix.m, chebfun, California.mat
%% ======================================================================
clear all;
rng(0);          % fixed seed so the random-matrix rows are reproducible

% Each row: {label, M, niter}
cases = {
    'Random 500',  10,  100
    'Random 1000', 10,  100
    'Poisson',     200, 300
    'California',  30,  100
};

nrows = size(cases,1);
results = cell(nrows,4);   % label, size, iters, err

for r = 1:nrows
    label = cases{r,1};  M = cases{r,2};  niter = cases{r,3};
    tol = 1e-14;

    % --- build A, b -----------------------------------------------------
    switch label
        case 'Random 500'
            N = 500;  A = rand(N); A = A/norm(A);
            b = rand(N,1); b = b/norm(b);
        case 'Random 1000'
            N = 1000; A = rand(N); A = A/norm(A);
            b = rand(N,1); b = b/norm(b);
        case 'Poisson'
            [A,b] = poisson_example(50); b = b/norm(b); N = length(A);
        case 'California'
            S = load('California.mat'); A = S.Problem.A; N = length(A);
            b = ones(N,1); b = b/norm(b);
    end

    exact = expm(full(A))*b;

    % --- Heaviside coefficient matrix ----------------------------------
    heav = @(t) 1 + 0*t;
    T = genCoeffMatrix(heav,M); T = sparse(T);
    w = [];
    for i=1:M, w(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2); end
    T(end-2:end,:) = sparse(3,M);

    % --- GMRES on the full star-linear system --------------------------
    Atimes = @(v) v - reshape((T/2*reshape(v,M,N))*A.',M*N,1);
    [x,flag,relres,iter,resvec] = ...
        gmres(Atimes,kron(b,w),[],tol,niter,[],[],kron(b,w));
    Star_Sol = sqrt(2)*x(1:M:end);

    err = norm(exact-Star_Sol)/norm(exact);

    results(r,:) = {label, N, iter(2), err};
end

%% print the table
fprintf('\n%-14s %6s %8s %18s\n','Matrix','Size','Iters','Err second star');
fprintf('%s\n', repmat('-',1,48));
for r = 1:nrows
    fprintf('%-14s %6d %8d %18.4e\n', results{r,:});
end
fprintf('\n');

%% ----------------------------------------------------------------------
function [L,v] = poisson_example(J)
% TREFETHEN, WEIDEMAN, and SCHMELZER 2006
h = 2/J; s = (-1+h:h:1-h)';
[xx,yy] = meshgrid(s,s);
x = xx(:); y = yy(:);
L = -gallery('poisson',J-1)/h^2;
I = speye((J-1)^2);
N = 32; theta = pi*(1:2:N-1)/N;
z = N*(.1309-.1194*theta.^2+.2500i*theta);
w = N*(-.1194*2*theta+.2500i);
c = (1i/N)*exp(z).*w;
u = (1-x.^2).*(1-y.^2).*exp(x);
v = zeros(size(u));
for k = 1:N/2
    v = v - c(k)*((z(k)*I-0.02*L)\u);
end
end
