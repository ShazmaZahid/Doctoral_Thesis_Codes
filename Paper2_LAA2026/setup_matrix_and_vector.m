function [Aop, b] = setup_matrix_and_vector(example_id, n, tmax)
%SETUP_MATRIX_AND_VECTOR  Build the test matrix and vector for each example.
%   example_id selects the benchmark problem; n sets its size; tmax scales A.
%   The returned vector b is normalized.

switch example_id
    case 1 % 2D Poisson matrix
        Aop = -gallery('poisson', n)*tmax;
        b = linspace(-1, 1, n^2)';

    case 2 % Complex tridiagonal matrix
        main_diag = 2i * ones(n, 1);
        off_diag  = -1i * ones(n-1, 1);
        A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
        Aop = zeros(n + 2);
        Aop(2:end-1, 2:end-1) = A;
        eps_val = 1e-13;
        Aop(1, 2) = eps_val; Aop(2, 1) = eps_val;
        Aop(end, end-1) = eps_val; Aop(end-1, end) = eps_val;
        Aop(1, end-1) = eps_val; Aop(end-1, 1) = eps_val;
        Aop = Aop * tmax;
        b = zeros(n + 2, 1); b(1) = 1;
        Aop = sparse(Aop);

    case 3 % Dense matrix with exponentially decaying eigenvalues
        lambda = exp(-linspace(0, 5, n));
        A = diag(lambda);
        Aop = A * tmax;
        b = rand(n, 1);
        [Q,~] = qr(rand(n,n));
        Aop = Q*sparse(Aop)*Q';

    case 4 % Tridiagonal symmetric Toeplitz
        c = [2; -1; zeros(n-2, 1)];
        r = [2, -1, zeros(1, n-2)];
        A = toeplitz(c, r);
        Aop = A * tmax;
        b = rand(n, 1);

    case 5 % Pentadiagonal nonsymmetric Toeplitz
        A = gallery('toeppen', n);
        Aop = A * tmax;
        b = rand(n, 1);

    case 6 % Dense matrix with eigenvalues at Chebyshev nodes
        lambda = cos((2*(1:n) - 1) * pi / (2*n));
        A = diag(lambda);
        Aop = A * tmax;
        b = rand(n, 1);
        [Q,~] = qr(rand(n,n));
        Aop = Q*sparse(Aop)*Q';

    case 7 % Schrodinger equation, FEM semidiscretization
        model = createpde();
        gd = [3,4,-2,2,2,-2,-2,-2,2,2]';
        ns = char('R1'); ns = ns';
        sf = 'R1';
        g = decsg(gd,sf,ns);
        geometryFromEdges(model,g);
        generateMesh(model,'Hmax',3/n);
        specifyCoefficients(model,'m',0,'d',1i,'c',1/2,'a',@V,'f',0);
        setInitialConditions(model,@u0);
        applyBoundaryCondition(model,'dirichlet', ...
            'Edge',1:model.Geometry.NumEdges,'u',0);
        xnode = model.Mesh.Nodes(1,:)';
        ynode = model.Mesh.Nodes(2,:)';
        u0vec = u0(struct('x',xnode,'y',ynode));
        FEM = assembleFEMatrices(model,"stiff-spring");
        Aop = -FEM.M\FEM.Ks*tmax;
        b = u0vec;
end

% Normalize the vector
b = b / norm(b);

end

%% Auxiliary functions for the Schrodinger example -----------------------
function V = V(x,~)
%V  Piecewise-constant potential for the Schrodinger equation.
n1 = 1; nr = numel(x.x);
V = zeros(n1,nr);
for i = 1:nr
    if x.x(i) >= -1 && x.x(i) <= 1 && x.y(i) >= -1 && x.y(i) <= 1
        V(1,i) = 0;
    else
        V(1,i) = 10.0;
    end
end
end

function u0 = u0(x)
%u0  Initial condition for the Schrodinger equation.
u0 = exp(-((x.x).^2 + (x.y).^2)/2);
end
