function [Star_Sol1] = run_star_method(Aop, b, M, T, w)
%RUN_STAR_METHOD  Base star-method (direct solve of the full system).
%   Solves the full star-linear system (I - A(x)T/2) x = b(x)w by backslash
%   and returns the endpoint solution.  For small matrices only.

    AT = kron(Aop, T);
    L  = kron(b, w);
    RAv = (speye(M * size(Aop,1)) - AT/2) \ L;   % solve the full system
    e1 = zeros(M,1); e1(1) = 1;
    K  = e1' * sqrt(2);
    Star_Sol1 = kron(speye(size(Aop,1)), K) * RAv;
end
