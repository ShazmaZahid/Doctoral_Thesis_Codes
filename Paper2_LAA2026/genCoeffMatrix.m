function [F, k] = genCoeffMatrix(f,M)
%GENCOEFFMATRIX Computes the coefficient matrix of f(t)*Theta(t-s), with a
%smooth function f and the Heaviside theta function Theta, in a basis of
%orthonormal Legendre polynomials
%INPUT:
%   f = function of interest
%   M = size of basis
%OUTPUT:
%   F = MxM matrix containing the Fourier coefficients of f(t)*Theta(t-s)
%NOTE:
%   Requires chebfun for the computation of the Legendre coefficients of f

c = cheb2leg(chebcoeffs(chebfun(f)));
k = length(c);
F = zeros(M);
for d = 0:length(c)-1
    Bd = genBasisMatrix(d,M);
    F = F + c(d+1)*Bd*sqrt(2)/sqrt(2*d+1);
end
end


function Bd = genBasisMatrix(d,M)
%GENBASISMATRIX generates the basis matrix of orthonormal Legendre
%polynomials times theta in a basis spanned by orthonormal Legendre
%polynomials
%INPUT:
%   d = degree of the orthonormal Legendre polynomial to be expanded
%   M = size of the basis
%OUTPUT:
%   Bd = MxM matrix containing the coefficients in the basis

Z = diag(ones(M-1,1),-1) + diag(-ones(M-1,1),1); Z(1,1) = 1; Z(M+1,M) = 1;
fact = sqrt(2*(0:M-1)+1)'./sqrt(2*(0:M-1)+1);   % Constant factor
[Firstcol,LastRow] = Hankelpart(d,M+1);
ToeplRow = Toeplitzpart(d,M+1);
H = hankel(Firstcol,LastRow); H = H(1:M,1:M+1);
T = toeplitz(ToeplRow);       T = T(1:M,1:M+1);

Bd = sqrt(2*d+1)/sqrt(8)*fact.*((H.*T)*Z);
end


function Firstrow = Toeplitzpart(d,m)
%TOEPLITZPART Generates Toeplitz matrix appearing in formula of Lemma 0.2
%INPUT:
%   d = degree of the orthonormal Legendre polynomial to be expanded
%   m = size of the basis
%OUTPUT:
%   Firstrow = vector of length m, the first row of the Toeplitz matrix

Firstrow = zeros(1,m);
for alpha = 0:m-1
    if mod(alpha+d,2) == 0 && alpha <= d      % d+alpha even
        temp = 1;
        count = 1;
        for j = 1:d+alpha
            if j <= (d-alpha)/2
                temp = temp*(((d+alpha)/2+j)/(2*j))^2;
                count = count + 1;
            end
            if count <= d
                temp = temp/2^2;
                count = count + 1;
            end
            if j <= alpha
                jj = alpha-j+1;               % reversing order for stability
                temp = temp*(d+jj)/(d-alpha+jj);
            end
        end
        Firstrow(alpha+1) = temp*2;
    end
end
end


function [Firstcol,LastRow] = Hankelpart(d,m)
%HANKELPART Generates Hankel matrix appearing in formula of Lemma 0.2
%INPUT:
%   d = degree of the orthonormal Legendre polynomial to be expanded
%   m = size of the basis
%OUTPUT:
%   Firstcol = vector of length m, the first column of the Hankel matrix
%   LastRow  = vector of length m, the last row of the Hankel matrix

Firstcol = zeros(m,1);
LastRow  = zeros(1,m);

for sumCol1 = 0:m-1
    if sumCol1 >= d && mod(d+sumCol1,2) == 0
        temp = 1/(d+sumCol1+1);
        for j = 1:d
            temp = temp*(-d+sumCol1+2*j)/(-d+sumCol1+2*j-1);
        end
        Firstcol(sumCol1+1) = temp;
    end
end

for sumRowm = m-1:2*m-2
    if mod(d+sumRowm,2) == 0
        temp = 1/(d+sumRowm+1);
        for j = 1:d
            temp = temp*(-d+sumRowm+2*j)/(-d+sumRowm+2*j-1);
        end
        LastRow(sumRowm-m+2) = temp;
    end
end
end
