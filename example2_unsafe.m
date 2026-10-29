%% Example 2, controller designed WITHOUT the safety constraint
% Same system and data as example2_safe.m, but Theorems 1 and 2 only
% enforce stability; the safety requirement h(x) = 1 - x1 >= 0 is ignored.
% The online step uses the least-squares solution of N*theta = e (Vert = N\e) instead
% of the vertices of the consistency set.

close all;clear;clc
yalmip('clear')

%% Settings
eps     = 1e-2;         % noise bound on the data
cc      = 1e-1;         % decrease rate in Theorem 1: V(x+) <= V(x) - cc*||x||^2
x0      = [-1; 1; 1];   % initial state
n_steps = 10;           % closed-loop steps

rng(3)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

%% System and data

% yalmip setup
sdpvar x1 x2 x3 x4 x5 x6
var = [x1;x2;x3];           % x
var_lift = [x4;x5;x6];      % \tilde{x} = x+: x4 = next(x1), x5 = next(x2), x6 = next(x3)
var_aug = [var; var_lift];  % [x; \tilde{x}]

% true system x+ = f(x) + g*u (see utils/getSystem.m)
eg = 4;
[f,g] = getSystem(eg,var);

% problem dimensions
Cons = setCons(eg);

% generate data, with input u \in [-1,1]
X = getDataRND(eg,eps,Cons.T);

% consistency set: data matrices
[A,xi] = getConsEx2(X,Cons);

% check enough data
N = [A; -A];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];

% keep nontrivial constraints
[N, e] = nontrivial_constraints(N, e);

% a single model: the least-squares solution of N*theta = e
Vert = N\e;

df  = Cons.df;

%% Theorem 1 without the safety constraint

dp = 2;     % degree of vecter polynomial p = [p1; p2; p3]
dv = 2;     % degree of lyapunov function v

% g = [2; -1; 1] hard-coded
[p1,c1] =   polynomial(var_aug, dp, 0);
[p2,c2] =   polynomial(var_aug, dp, 0);
[p3,c3] =   polynomial(var_aug, dp, 0);
[p3u,c3u] = polynomial(var_aug, dp, 0); % p3u = p3*u

% lyapunov function v
[v, cv] =  polynomial(var, dv, 1);          % 1 -> without constant
[vt,cvt] = polynomial(var_lift, dv, 1);

% multipler s (not used without the safety constraint)
[s,cs] = polynomial(var_aug, dp, 0);

% monomials of f: [x3, x2^2, x1*x2^2]
[~,~,vf] = polynomial(var,df,1);
vf = vf([3 6 12]);

% auxiliary varible Y
Y = zeros(1, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(i) = polynomial(var_aug, dv, 0);
end

% decrease of V for all f in the consistency set
r = [- kron([p1, p2, p3], vf')];
d = [v  - vt - cc*(var'*var) - [p1, p2, p3]*var_lift + p3u];

F = [
     coefficients(r - Y*N, var_aug) == 0,...
     (sos(-Y*e + d)):'tagged',...
     sos(Y),...
     sos(v - var'*var),...
     cv == cvt,...
     ];

opts = sdpsettings('solver','mosek','verbose',0);

[sol,sol_u,sol_Q] = solvesos(F, [], opts, [c1;c2;c3;c3u;cv;cvt;cs;coefficients(Y,var_aug)]);
sol

%% Extract the solution of Theorem 1
syms z1 z2 z3 z4 z5 z6
var = [z1;z2;z3];
var_aug = [z1;z2;z3;z4;z5;z6];

cv = value(cv);
vv = monomials(var, 1:dv);
V = cv'*vv;
V = matlabFunction(V);

VF = monomials(var,1:df);
VF = VF([3 6 12]);
VF = matlabFunction(VF);

% get system in syms
[f,g] = getSystem(eg,var);
FF = matlabFunction(f);

%% Theorem 2 without the safety constraint: compute the input online

x = x0;
UU = [];
XX = x;
VV = V(x(1),x(2),x(3));

for j = 1:n_steps
    fprintf('%d-th iteration \n',j)
    % yalmip setup
    yalmip('clear')
    sdpvar u

    % monomials of f at the current state
    vf = VF(x(1),x(2),x(3));

    % decrease V for this single model
    F = [];
    for i = 1:length(Vert)
        xt = kron(eye(3),vf')*Vert+g*u;
        F = [ F;
             V(xt(1),xt(2),xt(3)) <= V(x(1),x(2),x(3))-cc*(x'*x);
             ];
    end
    obj = u*u;
    opts = sdpsettings('solver','mosek','verbose', 0);
    optimize(F,obj,opts);

    % apply u to the true system
    u = value(u);
    xnext = FF(x(1),x(2),x(3))+g*u;
    dV = V(x(1),x(2),x(3)) - V(xnext(1),xnext(2),xnext(3));

    x = xnext;
    UU = [UU, u];
    XX = [XX, xnext];
    VV = [VV, V(xnext(1),xnext(2),xnext(3))];

end

plot3(XX(1,:),XX(2,:),XX(3,:))
hold on

% unsafe half-space boundary x1 = 1
unsafe1 = [1  1  1  1];
unsafe2 = [-2 2  2 -2];
unsafe3 = [-2 -2 2  2];
fill3(unsafe1,unsafe2,unsafe3,'r');

xlabel x1
ylabel x2
zlabel x3
xlim([-1;2])
caz =-15;
cel = 54;
view(caz,cel)
hold off

figure()
plot(1:length(VV),VV)
xlabel('k')
ylabel('V(x_k)')
legend('lyap')
