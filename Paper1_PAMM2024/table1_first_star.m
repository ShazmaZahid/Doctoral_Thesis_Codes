%% ======================================================================
%  TABLE 1  (PAMM 2024)  --  FIRST STAR-METHOD
%  Runs every case in one go and prints the table:
%     classical Krylov (expm)  vs  first star-method.
%
%  Requires on the path: Arnoldimgs.m, genCoeffMatrix.m, chebfun,
%                        California.mat
%% ======================================================================
clear all;
rng(0);          % fixed seed so the random-matrix rows are reproducible

% Each row: {label, N, M, k}.  A and b are built per case below.
cases = {
    'Random 500',  500,  10,  6
    'Random 500',  500,  10,  8
    'Random 500',  500,  10, 10
    'Random 1000', 1000, 10,  4
    'Random 1000', 1000, 10,  6
    'Random 1000', 1000, 10,  8
    'Poisson',     2401, 200, 100
    'Poisson',     2401, 200, 125
    'Poisson',     2401, 200, 150
    'California',  9664, 30,  10
    'California',  9664, 30,  15
    'California',  9664, 30,  22
};

nrows = size(cases,1);
results = cell(nrows,5);   % label, size, k, err_pade, err_star

for r = 1:nrows
    label = cases{r,1};  M = cases{r,3};  k = cases{r,4};

    % --- build A, b for this case ---------------------------------------
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

    % --- Arnoldi --------------------------------------------------------
    [Vk,Hk] = Arnoldimgs(A,b,k);
    Hk = Hk(1:k-1,1:k-1);

    % --- Pade Krylov baseline ------------------------------------------
    e1 = zeros(k-1,1); e1(1) = 1;
    Pade_Sol = Vk(:,1:k-1)*expm(Hk)*e1;

    % --- first star-method ---------------------------------------------
    heav = @(t) 1 + 0*t;
    T = genCoeffMatrix(heav,M); T = sparse(T);
    w = [];
    for i=1:M, w(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2); end
    T(end-2:end,:) = sparse(3,M);
    [U,S] = schur(full(T/2),'complex');
    Uw = (U'*w);
    Y = dlyap(S,Hk.', Uw*e1');
    X = U*Y;
    Star_Sol = sqrt(2)*Vk(:,1:k-1)*X(1,:).';

    err_pade = norm(exact-Pade_Sol)/norm(exact);
    err_star = norm(exact-Star_Sol)/norm(exact);

    results(r,:) = {label, N, k, err_pade, err_star};
end

%% print the table
fprintf('\n%-14s %6s %6s %18s %18s\n', ...
        'Matrix','Size','Iters','Err classical','Err first star');
fprintf('%s\n', repmat('-',1,66));
for r = 1:nrows
    fprintf('%-14s %6d %6d %18.4e %18.4e\n', results{r,:});
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
