function [A,xi] = getCons1(X,Cons) 
% 
% Output: 
%   A,B in (11)
%
% Input:
%   X: data matrix, columns represent [u,t,x1,x2,dx1,dx2]
%   df: degree of function f(x)
%   dg: degree of function g(x)
%   n: # of states
%   T: # of samples

kkk = 1;

df = Cons.df;
dg = Cons.dg;
n = Cons.n;
T = Cons.T;

lf = nchoosek(df+n, n) -kkk;	% length of monomial vector for f(x)
lg = nchoosek(dg+n, n);         % ... for g(x)

syms x1
vars = x1;

p = monomials(vars, kkk:df); 	% phi
g = monomials(vars, 0:dg);      % gamma

p_fun = matlabFunction(p);
g_fun = matlabFunction(g);

A = []; B = [];
p_data  = zeros(T, lf);
ug_data = zeros(T, lg);

ind = randperm(size(X,1),T);
X = X(ind,:);

for i = 1:T
    p_data(i,:) = feval(p_fun, X(i,3));
    
    if dg == 0
        ug_data(i,:) = X(i, 1);
    elseif dg > 0
        ug_data(i, :) = X(i, 1)*feval(g_fun, X(i, 3));
    end
    
    A = [A; kron(eye(n),p_data(i,:))];
    B = [B; kron(eye(n),ug_data(i,:))];
end

dx = X(1:T,5)'-1*X(1:T,1)';
xi = dx(:);

end