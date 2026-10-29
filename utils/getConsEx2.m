function [A,xi] = getConsEx2(X,Cons)
% GETCONSEX2  Data matrices of the consistency set for Example 2 (3D).
%   [A, xi] = getConsEx2(X, Cons) returns A and xi such that every f
%   compatible with the data satisfies |A*theta - xi| <= eps, where theta
%   are the coefficients of f in the dictionary [x3, x2^2, x1*x2^2] and
%   xi = x+ - g*u with the known g = [2; -1; 1].
%
%   X     data matrix, one sample per row: [u, t, x1, x2, x3, x1+, x2+, x3+]
%   Cons  problem dimensions (n, T, df)
%
%   Note: the rows of X are shuffled with randperm (uses the global random
%   stream).

dmin = 1;   % lowest monomial degree in f (f(0) = 0)

df = Cons.df;
n = Cons.n;
T = Cons.T;

syms x1 x2 x3
vars = [x1,x2,x3];

p = monomials(vars, dmin:df); 	% phi
p = p([3 6 12]);                % [x3, x2^2, x1*x2^2]
lf = length(p);

p_fun = matlabFunction(p);

A = [];
p_data  = zeros(T, lf);

ind = randperm(size(X,1),T);
X = X(ind,:);

for i = 1:T
    p_data(i,:) = feval(p_fun, X(i,3), X(i,4), X(i,5));
    A = [A; kron(eye(n),p_data(i,:))];
end

dx = X(1:T,[6 7 8])' - [2;-1;1].*X(1:T,1)'; % g = [2; -1; 1] is known
xi = dx(:);

end
