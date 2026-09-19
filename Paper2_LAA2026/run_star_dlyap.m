function [Star_Sol] = run_star_dlyap(Aop, b, S, U, Uw)
%RUN_STAR_DLYAP  Schur-based star-method (direct variant).
%   Solves the transformed Stein equation with dlyap, using the precomputed
%   Schur factor S of T_M/2 (so dlyap skips its internal Schur step), and
%   returns the endpoint solution.
%       S, U : Schur factors of T_M/2  (T_M/2 = U S U^H)
%       Uw   : U' * w, with w the Legendre weight vector phi_M(-1)

    Y = dlyap(S, Aop.', Uw*b.');   % solve the (triangular) Stein equation
    X = U * Y;                     % transform back
    Star_Sol = sqrt(2) * X(1, :).';
end
