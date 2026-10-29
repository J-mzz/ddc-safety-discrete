%% Example 2 of "Data-Driven Safe Control of Discrete-Time Non-Linear Systems"
% Safe controller for a 3D polynomial system with a convex (half-space)
% safety requirement (Fig. 2).
%
% True system (used to generate the data and to simulate):
%   x+ = [x2^2 + x3; x1*x2^2 + x3; 0] + [2; -1; 1]*u
% Safe set: h(x) = 1 - x1 >= 0
% Prior:    f is spanned by [x3, x2^2, x1*x2^2], g = [2; -1; 1] is known;
%           8 samples with ||w||_inf <= 0.01
%
% Steps
%   1. Build the consistency set of f from noisy data and the bounding box
%      of its vertices.
%   2. Theorem 1: SOS program for a Lyapunov function V that certifies safe
%      stabilization of every system in the consistency set.
%   3. Theorem 2: at every step, compute the input online by a convex
%      program that decreases V and keeps x+ safe for all box vertices.
%
% Figures: 3D trajectory with the unsafe half-space, and V along the
% trajectory. example2_unsafe.m designs the controller without the safety
% constraint.

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

% find vertices
Vert = con2vert(N,e);
oldVert = Vert;

% bounding box for vertices
vmax = max(Vert,[],1);
vmin = min(Vert,[],1);
newVert = [];
vertlen = length(vmax);
for number_vert = 1:2^vertlen
    binStr = dec2bin(number_vert,vertlen);
    vert_new = zeros(1,vertlen);
    for jj = 1:vertlen
        if strcmp(binStr(jj),'1')
            vert_new(jj) = vmax(jj);
        else
            vert_new(jj) = vmin(jj);
        end
    end
    newVert = [newVert;vert_new];
end
Vert = newVert;

%% Safe set h(x) = 1 - x1 >= 0 (half-space)
h  = 1-x1;
ht = 1-x3;

df  = Cons.df;

%% Theorem 1: SOS program for the Lyapunov function V

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

% multipler s: sos( H(xt)-s1*H(x) )
[s,cs] = polynomial(var_aug, dp, 0);

% monomials of f: [x3, x2^2, x1*x2^2]
[~,~,vf] = polynomial(var,df,1);
vf = vf([3 6 12]);

% auxiliary varible Y
Y = zeros(2, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(1,i) = polynomial(var_aug, dv, 0);
    Y(2,i) = polynomial(var_aug, dv, 0);
end

% decrease of V (first row) and safety (second row) for all f in the
% consistency set
r = [- kron([p1, p2, p3], vf');
     - kron([p1, p2, p3], vf')];
d = [v  - vt - cc*(var'*var) - [p1, p2, p3]*var_lift + p3u;
     ht - s*h                - [p1, p2, p3]*var_lift + p3u];

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

H = 1-z1;
H = matlabFunction(H);

VF = monomials(var,1:df);
VF = VF([3 6 12]);
VF = matlabFunction(VF);

% get system in syms
[f,g] = getSystem(eg,var);
FF = matlabFunction(f);

%% Theorem 2: compute the input online

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

    % decrease V and stay safe for every vertex of the bounding box
    F = [];
    for i = 1:length(Vert)
        xt = kron(eye(3),vf')*Vert(i,:)'+g*u;
        F = [ F;
             V(xt(1),xt(2),xt(3)) <= V(x(1),x(2),x(3))-cc*(x'*x);
             H(xt(1)) >= 0;
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
