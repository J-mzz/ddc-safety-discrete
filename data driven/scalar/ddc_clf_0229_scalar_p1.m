close all;clear;clc
yalmip('clear')

rng(1)

%% system and data generation  system with dim 1

% yalmip setup
sdpvar x1 x3
var = x1;              % x
var_lift = x3;         % \tilde{x}: x3 = next(x1), x4 = next(x2)
var_aug = [var; var_lift];  % [x; \tilde{x}]

% define DT system: x+ = f(x) + g(x)*u(x) + noise
eg = 3;                     % 1 (linear) or 2 (nonlinear) or 3 (scalar)
[f,g] = getSystem(eg,var);  % f = -2x; g = 1

% setup
Cons = setCons(eg);

% noise level
eps = 1e-1;                    % noise level for robust data driven
epsw = 0;                   % noise level of disturbance

% generate data, with input u \in [-1,1]
X = getDataRND(eg,eps,Cons.T);

% consistency set
[A,xi] = getCons1(X,Cons);    % for eg = 3

% check enough data
N = [A; -A];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];

% keep nontrivial constraints
[N, e] = nontrivial_constraints(N, e);

%% try to solve
n   = Cons.n;
T   = Cons.T;
df  = Cons.df;
dg  = Cons.dg;

%% Solve Thm 1

dp = 2;     % degree of vecter polynomial p = [p1; p2]
dv = 2;     % degree of lyapunov function v

% g = [0; 1] hard-coded
[p2,c2] =   polynomial(var_aug, dp, 0);
[p2u,c2u] = polynomial(var_aug, dp, 0); % p2u = p2*u

% lyapunov function v
[v, cv] =  polynomial(var, dv, 1);          % 1 -> without constant
[vt,cvt] = polynomial(var_lift, dv, 1);

% vp and vf
[~,~,vf] = polynomial(var,df,1);

% auxiliary varible Y
Y = zeros(1, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(i) = polynomial(var_aug, dv, 0);
end

% for stable and safety set P2
c3 = 1e-4;
r = - kron(p2, vf');
d = [v  - vt - c3*(var'*var) - p2*var_lift + p2u];

F = [ 
     coefficients(r - Y*N, var_aug) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     sos(v - var'*var),...
     cv == cvt,...
     ];

opts = sdpsettings('solver','mosek','verbose',0);

sol = solvesos(F, [], opts, [c2;c2u;cv;cvt;coefficients(Y,var_aug)])

%% extract solution Thm1
syms z1 z3
var = z1;
var_aug = [z1;z3];

cv = value(cv);   %  cv = cv.*(cv>=1e-6 | cv<=-1e-6);    % clean small values
vv = monomials(var, 1:dv);
V = cv'*vv;
V = matlabFunction(V);

c2  = value(c2);  %  c2  = c2.*(c2>=1e-6 | c2<=-1e-6);
c2u = value(c2u); %  c2u = c2u.*(c2u>=1e-6 | c2u<=-1e-6);
c2u(1) = 0;

vp = monomials(var_aug, 0:dp);
P2 = c2'*vp;
P2U = c2u'*vp;

P2 = matlabFunction(P2);
P2U = matlabFunction(P2U);

VF = monomials(var,1:df);
VF = matlabFunction(VF);

VG = 1;

% get system in syms
[f,g] = getSystem(eg,var); 
FF = matlabFunction(f);

% check if V is sos
figure()
fc = fcontour(V,[-2 2 -2 2]);
fc.LevelList = [0:20];
colorbar
axis square
hold on

%% Jared: using V, rho, find min-norm u(x),

%this should be an sos problem, not a moment problem
%but sos feasibility might not be guaranteed at finite degree?
% initial point 
x = 1;

% yalmip setup
yalmip('clear')

sdpvar x3
var_lift = x3;
u = sdpvar(1,1);

Y = zeros(1, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(i) = polynomial(var_lift, dv+2, 0);
end

% multiplier p
p2  = P2(x, x3);

% vf
vf = VF(x);

% constraints:
r = -p2*vf';
d = V(x)-V(x3)-c3*x^2 - p2*(x3-u);

F = [ 
     coefficients(r - Y*N, var_lift) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     ];

obj = u^2;

sol = solvesos(F, obj, opts, [u; coefficients(Y,var_lift)])

u_rec = value(u);

% P = msdp(min(obj),F,2);

% mset('yalmip',true)
% mset(sdpsettings('solver','mosek'))

% [status,obj] = msol(P)

%% check 
% x3 = double(x3)
% u = double(u)
% x_curr = x
% x_next = FF(x)+u
% dV = V(x) - V(x_next)
