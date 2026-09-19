%% ======================================================================
%  SCHUR-BASED STAR-METHODS               Pozza & Zahid (2026)
%  Numerical experiments of Paper 2 (Chapter 4 of the thesis, Section 4.3):
%  the star-method and its Arnoldi variant, compared with expv and
%  expmv_tspan, over ten benchmark problems.  Reproduces Tables 4.3-4.5.
%
%  Requires on the MATLAB path:
%     setup_matrix_and_vector.m, Arnoldimgs.m, genCoeffMatrix.m,
%     run_star_method.m, run_star_dlyap.m, run_star_krylov.m,
%     expv.m, expmv_tspan.m  (third party, see README),  and chebfun.
%
%  Select the experiment with  test_id  below (1-10).  Timings and errors
%  are written to a dated text file.
%% ======================================================================
clear all;

tmax  = 1;
t0    = 0;
nsample = 100;   % repetitions for average timing (expm excluded)

dat = date;
fileID = fopen(['data_comp_time',datestr(dat),'.txt'],'w');
fprintf(fileID,['Computational time - ', datestr(dat),'\n']);
fprintf(fileID,'nsample= %d \n \n',nsample);

for test_id = 5     % choose 1..10

    switch test_id
        case 1   % 2D Poisson matrix
            example_id = 1; n = 50;
            k = 35; M = 22; q = M; star_original = 0; star_dlyap = 1; tmax = 4;
        case 2   % complex tridiagonal matrix
            example_id = 2; n = 100;
            k = 17; M = 7;  q = M; star_original = 1; star_dlyap = 1; tmax = 8;
        case 3   % dense matrix, decaying eigenvalues (n = 2000)
            example_id = 3; n = 2000;
            k = 17; M = 13; q = M; star_original = 0; star_dlyap = 0; tmax = 4;
        case 4   % dense matrix, decaying eigenvalues (n = 20)
            example_id = 3; n = 20;
            k = 19; M = 15; q = M; star_original = 1; star_dlyap = 1; tmax = 4;
        case 5   % tridiagonal Toeplitz matrix
            example_id = 4; n = 100;
            k = 30; M = 25; q = M; star_original = 1; star_dlyap = 1; tmax = 4;
        case 6   % pentadiagonal Toeplitz matrix
            example_id = 5; n = 1000;
            k = 80; M = 38; q = M; star_original = 1; star_dlyap = 1; tmax = 2;
        case 7   % dense matrix, eigenvalues at Chebyshev nodes
            example_id = 6; n = 500;
            k = 20; M = 12; q = M; star_original = 1; star_dlyap = 1; tmax = 4;
        case 8   % Schrodinger FEM semidiscretization (size 841)
            example_id = 7; n = 10;
            k = 70; M = 50;  q = M; star_original = 0; star_dlyap = 1; tmax = 1/10;
        case 9   % Schrodinger FEM semidiscretization (size 1917)
            example_id = 7; n = 15;
            k = 120; M = 150; q = M; star_original = 0; star_dlyap = 1; tmax = 1/10;
        case 10  % Schrodinger FEM semidiscretization (size 3437)
            example_id = 7; n = 20;
            k = 180; M = 190; q = M; star_original = 0; star_dlyap = 1; tmax = 1/10;
    end

    %% Setup matrix and vector
    [Aop, b] = setup_matrix_and_vector(example_id, n, 1);

    fprintf(fileID,'example_id = %d, n = %d , size = %d , k = %d , M = %d , tmax = %d \n \n', ...
            example_id,n,length(Aop),k,M,tmax);

    %% Reference solution
    tic;
    exact = expm(Aop*tmax) * b;
    expm_time = toc;
    fprintf(fileID,['Exact expm computation time: ', num2str(expm_time),'\n']);

    %% expmv_tspan (applied to A)
    tic;
    for j=1:nsample
        [X_tspan, tvals, mv] = expmv_tspan(Aop, b, t0, tmax, q);
    end
    expmv_tspan_time = toc/nsample;
    fprintf(fileID,['expmv_TSPAN approximation time: ', num2str(expmv_tspan_time),'\n']);
    Expmv_Sol = X_tspan(:, end);
    error_method_expmv_tspan = norm(Expmv_Sol - exact) / norm(exact);
    fprintf(fileID,['Relative error (expmv_TSPAN) ', num2str(error_method_expmv_tspan),'\n']);

    %% expv (applied to A)
    t = tmax;
    tic
    for j=1:nsample
        [w_expv, err] = expv(t, Aop, b);
    end
    expv_time = toc/nsample;
    fprintf(fileID,['expv approximation time: ', num2str(expv_time),'\n']);
    error_expmv = norm(w_expv - exact) / norm(exact);
    fprintf(fileID,['Relative error (w_expv): ', num2str(error_expmv ),'\n']);

    %% Arnoldi projection (shared)
    tic;
    for j=1:nsample
        [Vk, Hk] = Arnoldimgs(Aop, b, k);
        Hk = Hk(1:k-1, 1:k-1);
    end
    arnoldi_time = toc/nsample;
    fprintf(fileID,['Arnoldi process time: ', num2str(arnoldi_time),'\n']);

    %% Heaviside coefficient matrix and precomputed Schur factor
    heav = @(t) 1 + 0*t;
    T = genCoeffMatrix(heav, M);
    T(end,:) = zeros(1,M);              % one-row truncation
    w = zeros(M, 1);
    for i = 1:M
        w(i) = (-1)^(i-1) * sqrt((2*(i-1)+1)/2);
    end
    [U, S] = schur(full(T/2), 'complex');
    Uw = U' * w;

    %% Star-method with Arnoldi
    tic
    for j=1:nsample
        [Star_Sol] = run_star_krylov(Vk, Hk, S*tmax, U, Uw);
    end
    star_krylov_time = toc/nsample;
    fprintf(fileID,['Star + Krylov approximation time: ', num2str(star_krylov_time+arnoldi_time),'\n']);
    error_method_star_krylov = norm(Star_Sol - exact) / norm(exact);
    fprintf(fileID,['Relative error (Star Krylov): ', num2str(error_method_star_krylov),'\n']);

    %% Classical Krylov + Pade (baseline)
    tt = linspace(t0,tmax,q); nb = norm(b);
    tic
    for j=1:nsample
        for i=1:q
            tmp = expm(Hk*tt(i));
            Pade_Sol = Vk(:,1:k-1)*(tmp(:,1)*nb);
        end
    end
    pade_krylov_time = toc/nsample;
    fprintf(fileID,['Pade + Krylov approximation time: ', num2str(pade_krylov_time+arnoldi_time),'\n']);
    error_pade = norm(Pade_Sol - exact) / norm(exact);
    fprintf(fileID,['Relative error (Pade Krylov): ', num2str(error_pade),'\n']);

    %% expmv_tspan + Arnoldi (baseline)
    e1nb = eye(k-1,1)*nb;
    tic
    for j=1:nsample
        [X_tspan, tvals, mv] = expmv_tspan(Hk, e1nb, t0, tmax, q);
        XK_tspan = Vk(:,1:k-1)*X_tspan;
    end
    Expmv_Sol = XK_tspan(:, end);
    expmvt_krylov_time = toc/nsample;
    fprintf(fileID,['expmv_TSPAN + Krylov approximation time: ', num2str(expmvt_krylov_time+arnoldi_time),'\n']);
    error_expmvt_krylov = norm(Expmv_Sol - exact) / norm(exact);
    fprintf(fileID,['Relative error (expmv_TSPAN + Krylov) ', num2str(error_expmvt_krylov),'\n']);

    %% Star-method, direct (full system, backslash)
    if star_original
        tic
        for j=1:nsample
            [Star_Sol1] = run_star_method(Aop*tmax, b, M, T, w);
        end
        starproduct_time = toc/nsample;
        fprintf(fileID,['Star approximation time: ', num2str(starproduct_time),'\n']);
        error_method_star = norm(Star_Sol1 - exact) / norm(exact);
        fprintf(fileID,['Relative error (Star Solution): ', num2str(error_method_star),'\n']);
    end

    %% Star-method, direct via Schur + dlyap (full A)
    if star_dlyap
        tic
        for j=1:nsample
            [Star_Sol2] = run_star_dlyap(Aop*tmax, b, S, U, Uw);
        end
        starproductdlyap_time = toc/nsample;
        fprintf(fileID,['Star dlyap approximation time: ', num2str(starproductdlyap_time),'\n']);
        error_method_star_dlyap = norm(Star_Sol2 - exact) / norm(exact);
        fprintf(fileID,['Relative error (Star dlyap Solution): ', num2str(error_method_star_dlyap),'\n']);
    end

    fprintf(fileID,'\n ----------- \n');
end

fclose(fileID);
