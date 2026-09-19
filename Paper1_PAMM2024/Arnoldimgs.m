function [V,H] = Arnoldimgs(A,s,m)
% Modified Arnoldi algorithm.
% [V,H] = Arnoldimgs(A,s,m), for a square matrix A and vector s of
% corresponding size, it produces the orthogonal basis V for the Krylov
% subspace K_m(A,s) = span{s, As, ..., A^{m-1}s}, and the m x m-1 matrix H
% so that: A V(:,1:m-1) = V H.
% Note that H(1:m-1,1:m-1) is an upper-Hessenberg matrix.
% The algorithm is implemented using the modified Gram-Schmidt recurrences.

tol = 1e-15;                          % Breakdown tolerance
V(:,1) = s/norm(s);
for i = 1:m-1
    z = A*V(:,i);
    for j = 1:i
        H(j,i) = V(:,j)'*z;
        z = z - V(:,j) * H(j,i);      % Modified Gram-Schmidt orthogonalization
    end
    H(i+1,i) = norm(z);
    if norm(z) < tol                  % Breakdown check
        return
    else
        V(:,i+1) = z/H(i+1,i);        % Normalization
    end
end
