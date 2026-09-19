function [Star_Sol] = run_star_krylov(Vk, Hk, S, U, Uw)
%RUN_STAR_KRYLOV  Schur-based star-method with Arnoldi (reduced variant).
%   Projects onto the Krylov subspace (basis Vk, Hessenberg Hk) and solves
%   the reduced Stein equation with dlyap, using the precomputed Schur
%   factor S of T_M/2.  Returns the endpoint solution.
%       Vk, Hk : Arnoldi basis and reduced Hessenberg matrix
%       S, U   : Schur factors of T_M/2
%       Uw     : U' * w, with w the Legendre weight vector phi_M(-1)

    Y = dlyap(S, Hk.', [Uw, zeros(size(Uw,1), size(Hk,1)-1)]);
    X = U * Y;                     % transform back
    Star_Sol = sqrt(2) * Vk(:, 1:size(Hk,1)) * X(1, :).';
end
