function [A,xi] = getCons3(X,Cons) 
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

syms x1 x2 x3
vars = [x1,x2,x3];

p = monomials(vars, kkk:df); 	% phi
p = p([3 6 12]);
lf = length(p);
g = monomials(vars, 0:dg);      % gamma

p_fun = matlabFunction(p);

g_fun = matlabFunction(g);

A = []; B = [];
p_data  = zeros(T, lf);
ug_data = zeros(T, lg);

ind = randperm(size(X,1),T);
X = X(ind,:);

for i = 1:T
    p_data(i,:) = feval(p_fun, X(i,3), X(i,4), X(i,5));
    
    if dg == 0
        ug_data(i,:) = X(i, 1);
    elseif dg > 0
        ug_data(i, :) = X(i, 1)*feval(g_fun, X(i, 3), X(i, 4), X(i,5));
    end
    
    A = [A; kron(eye(n),p_data(i,:))];
    B = [B; kron(eye(n),ug_data(i,:))];
end

dx = X(1:T,[6 7 8])' - [2;-1;1].*X(1:T,1)'; % g=[0;1] fixed
xi = dx(:);

end