%% Example 1 of "Data-Driven Safe Control of Discrete-Time Non-Linear Systems"
% Safe controller for a linear system with a convex safety requirement
% (blue trajectory in Fig. 1).
%
% True system (used to generate the data and to simulate):
%   x+ = [x2; -2*x1 - 5*x2] + [0; 1]*u,   open-loop unstable
% Safe set: h(x) = 4 - x1^2 - (x2 - 1)^2 >= 0
% Prior:    f is linear, g = [0; 1] is known; 40 samples with ||w||_inf <= 0.1
%
% Steps
%   1. Build the consistency set of f from noisy data (polytope N*f <= e)
%      and its vertices.
%   2. Theorem 1: SOS program for a Lyapunov function V that certifies safe
%      stabilization of every system in the consistency set.
%   3. Theorem 2: at every step, compute the input online by a convex
%      program that decreases V and keeps x+ safe for all vertices.
%
% Figures: trajectory with the level sets of V and the safe-set boundary,
% and V along the trajectory. example1_unsafe.m designs the controller
% without the safety constraint (orange trajectory in Fig. 1).

close all;clear;clc
yalmip('clear')

%% Settings
eps     = 1e-1;     % noise bound on the data
c3      = 1e-1;     % decrease rate in Theorem 1: V(x+) <= V(x) - c3*||x||^2
x0      = [2; 1];   % initial state
n_steps = 40;       % closed-loop steps

rng(12)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

%% System and data

% yalmip setup
sdpvar x1 x2 x3 x4
var = [x1;x2];              % x
var_lift = [x3;x4];         % \tilde{x} = x+: x3 = next(x1), x4 = next(x2)
var_aug = [var; var_lift];  % [x; \tilde{x}]

% true system x+ = f(x) + g*u (see utils/getSystem.m)
eg = 1;
[f,g] = getSystem(eg,var);

% problem dimensions
Cons = setCons(eg);

% generate data, with input u \in [-1,1]
X = getDataRND(eg,eps,Cons.T);

% consistency set: data matrices
[A,xi] = getConsEx1(X,Cons);

% check enough data
N = [A; -A];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];

% keep nontrivial constraints
[N, e] = nontrivial_constraints(N, e);
Vert = con2vert(N,e);

%% Safe set h(x) >= 0: disk with squared radius r0 and center c0 (convex)
r0 = 4;         % squared radius
c0 = [0; 1];	% center

h  = -r0 + (x1-c0(1))^2 + (x2-c0(2))^2;
ht = -r0 + (x3-c0(1))^2 + (x4-c0(2))^2;
h = -h;         % h(x)  >= 0
ht = -ht;       % h(x+) >= 0

df  = Cons.df;

%% Theorem 1: SOS program for the Lyapunov function V

dp = 2;     % degree of vecter polynomial p = [p1; p2]
dv = 2;     % degree of lyapunov function v

% g = [0; 1] hard-coded
[p1,c1] =   polynomial(var_aug, dp, 0);
[p2,c2] =   polynomial(var_aug, dp, 0);
[p2u,c2u] = polynomial(var_aug, dp, 0); % p2u = p2*u

% lyapunov function v
[v, cv] =  polynomial(var, dv, 1);          % 1 -> without constant
[vt,cvt] = polynomial(var_lift, dv, 1);

% multipler s: sos( H(xt)-s1*H(x) )
[s,cs] = polynomial(var_aug, dp, 0);

% monomials of f
[~,~,vf] = polynomial(var,df,1);

% auxiliary varible Y
Y = zeros(2, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(1,i) = polynomial(var_aug, dv, 0);
    Y(2,i) = polynomial(var_aug, dv, 0);
end

% decrease of V (first row) and safety (second row) for all f in the
% consistency set
r = [- kron([p1, p2], vf');
     - kron([p1, p2], vf')];
d = [v  - vt - c3*(var'*var) - [p1, p2]*var_lift + p2u;
     ht - s*h                - [p1, p2]*var_lift + p2u];

F = [
     coefficients(r - Y*N, var_aug) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     sos(v - var'*var),...
     cv == cvt,...
     ];

opts = sdpsettings('solver','mosek','verbose',0);

sol = solvesos(F, [], opts, [c1;c2;c2u;cv;cvt;cs;coefficients(Y,var_aug)])

%% Extract the solution of Theorem 1
syms z1 z2 z3 z4
var = [z1;z2];
var_aug = [z1;z2;z3;z4];

cv = value(cv);
vv = monomials(var, 1:dv);
V = cv'*vv;
V = matlabFunction(V);

H = -r0 + (z1-c0(1))^2 + (z2-c0(2))^2;
H = -H;
H = matlabFunction(H);

VF = monomials(var,1:df);
VF = matlabFunction(VF);

% get system in syms
[f,g] = getSystem(eg,var);
FF = matlabFunction(f);

% level sets of V
figure()
fc = fcontour(V,[-2 2 -2 2]);
fc.LevelList = [0:20];
axis square
hold on

%% Theorem 2: compute the input online

x = x0;
UU = [];
XX = x;
VV = V(x(1),x(2));

for j = 1:n_steps
    fprintf('%d-th iteration \n',j)
    % yalmip setup
    yalmip('clear')
    sdpvar u

    % monomials of f at the current state
    vf = VF(x(1),x(2));

    % decrease V and stay safe for every vertex of the consistency set
    F = [];
    for i = 1:length(Vert)
        xt = kron(eye(2),vf')*Vert(i,:)'+g*u;
        F = [ F;
             V(xt(1),xt(2)) <= V(x(1),x(2))-c3*(x'*x);
             H(xt(1),xt(2)) >= 0;
             ];
    end
    obj = u*u;

    opts = sdpsettings('solver','mosek','verbose', 0);
    sol = optimize(F,obj,opts);
    if sol.problem
        error('Online problem at step %d: %s', j, sol.info)
    end

    % apply u to the true system
    u = value(u);
    xnext = FF(x(1),x(2))+g*u;
    dV = V(x(1),x(2)) - V(xnext(1),xnext(2));

    x = xnext;
    UU = [UU, u];
    XX = [XX, xnext];
    VV = [VV, V(xnext(1),xnext(2))];
end

plot(XX(1,:),XX(2,:))

% safe-set boundary h(x) = 0
theta = (0:60)/30*pi;
X_unsafe = c0(1) + sqrt(r0)*cos(theta);
Y_unsafe = c0(2) + sqrt(r0)*sin(theta);
plot(X_unsafe,Y_unsafe,'r')
ylim([-2 2])
xlim([-2 2])
axis square
hold off

figure()
plot(1:length(VV),VV)
legend('lyap')
