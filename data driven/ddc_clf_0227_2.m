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
[f,g] = getSystem(eg,var); 

% setup
Cons = setCons(eg);

% noise level
eps = 1e-1;                    % noise level for robust data driven
epsw = 0;                   % noise level of disturbance

% generate data, with input u \in [-1,1]
X = getDataRND(eg,eps,Cons.T);

% consistency set
if eg == 3
    [A,B,xi] = getCons2(X,Cons);    % for eg = 3
else
    [A,B,xi] = getCons(X,Cons);     % for eg = 1,2
end

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
syms z1
var = z1;

cv = value(cv);     cv = cv.*(cv>=1e-6 | cv<=-1e-6);    % clean small values
vv = monomials(var, 1:dv);
V = cv'*vv;
V = matlabFunction(V);

VF = monomials(var,1:df);
VF = matlabFunction(VF);

VG = 1;

% get system in syms
[f,g] = getSystem(eg,var); 
FF = matlabFunction(f);

% % check if V is sos
% figure()
% fc = fcontour(V,[-2 2 -2 2]);
% fc.LevelList = [0:20];
% colorbar
% axis square
% hold on

%% Thm2 step 1

% initial point 
x = 1;

% yalmip setup
yalmip('clear')
sdpvar alpha
sdpvar x3
var = x3;              % tilde x

% multiplier p
[p2,c2] =   polynomial(var, dp, 0);
[p2u,c2u] = polynomial(var, dp, 0);     % p2u = psi^T*g

% multipler s
[s,cs] = polynomial(var, dp, 0);

% vp and vf
vf = VF(x);
vg = VG;

% auxiliary varible Y
Y = zeros(1, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(i) = polynomial(var, dv, 0);
end

% for stable set P2
r = - kron(p2, vf');
d = alpha - V(x3) - p2*var + p2u - s;


F = [ 
     coefficients(r - Y*N, var) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     sos(s),...
     alpha >= V(x),...
     ];

obj = alpha;
opts = sdpsettings('solver','mosek','verbose',1);

sol = solvesos(F, obj, [], [c2;c2u;cv;cs;coefficients(Y,var);alpha])

alpha = value(alpha);
%}

%% extract solution thm2 step 1: fixing all except tilde x 

cY = value(coefficients(Y,var));
cY = cY.*(cY>=1e-6 | cY<=-1e-6);

syms z3
var = z3;

% clean small values
c2 = value(c2);     c2 = c2.*(c2>=1e-9 | c2<=-1e-9);
c2u = value(c2u);   c2u = c2u.*(c2u>=1e-6 | c2u<=-1e-6);
cs = value(cs);

v2 =  monomials(var, 0:dp);
v2u = monomials(var, 0:dp);
vs = monomials(var, 0:dp);

P2 = c2'*v2; 
P2U = c2u'*v2u;
S = cs'*vs;

vY = monomials(var, 0:dv);
nY = length(vY);
for i = 1:size(N,1)
    YY{i} = cY((1:nY)+nY*(i-1))'*vY;
    YY{i} = matlabFunction(YY{i});
end

P2 = matlabFunction(P2);
P2U = matlabFunction(P2U);
S = matlabFunction(S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixed everything except tilde x, with gloptipoly

% yalmip setup
yalmip('clear')

mpol x3
var = x3;              % \tilde x

% multiplier p
p2  = P2(x3);
p2u = P2U(x3);

% lyapunov function v
vt = V(x3);

% vp and vf
vf = VF(x);
vg = VG;

% new
Y = {};
for i = 1:size(N,1)
    Y{i} = YY{i}(x3);
end

% for stable set P2
r = - p2*vf';
d = V(x)  - V(x3) - c3*(x^2)- p2*var + p2u - S(x3);
% d = alpha  - V(x3) - p2*var + p2u - S(x3);

Ye = 0;
for i = 1:length(Y)
    Ye = Ye+Y{i}*e(i);
end

F = [ 
     -Ye + d >= 0,...
     ];

obj = V(x3);
P = msdp(max(obj),F,4);

mset('yalmip',true)
mset(sdpsettings('solver','mosek'))

[status,obj] = msol(P)

x3 = double(x3)

%% result

u = P2U(x3)/P2(x3)
xcurr = x
Vcurr = V(x)
alpha
xnext = FF(x) + g*u
Vnext = V(xnext)