close all;clear;clc
yalmip('clear')

rng(1)

%% system and data generation  system with dim 1

% yalmip setup
sdpvar x1 x2 x3 x4
var = [x1;x2];              % x
var_lift = [x3;x4];         % \tilde{x}: x3 = next(x1), x4 = next(x2)
var_aug = [var; var_lift];  % [x; \tilde{x}]

% define DT system: x+ = f(x) + g(x)*u(x) + noise
eg = 2;                     % 1 (linear) or 2 (nonlinear) or 3 (scalar)
[f,g] = getSystem(eg,var); 

% setup
Cons = setCons(eg);

% noise level
eps = 1e-2;                    % noise level for robust data driven
epsw = 0;                   % noise level of disturbance

% generate data, with input u \in [-1,1]
X = getDataRND(eg,eps,Cons.T);

% consistency set
[A,xi] = getCons2(X,Cons);    % for eg = 3

% check enough data
N = [A; -A];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];

% keep nontrivial constraints
[N, e] = nontrivial_constraints(N, e);
Vert = con2vert(N,e);
%% try to solve
n   = Cons.n;
T   = Cons.T;
df  = Cons.df;
dg  = Cons.dg;

%% Solve Thm 1

dp = 2;     % degree of vecter polynomial p = [p1; p2]
dv = 2;     % degree of lyapunov function v

% g = [0; 1] hard-coded
[p1,c1] =   polynomial(var_aug, dp, 0);
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
c3 = 1e-2;
r = - kron([p1, p2], vf');
d = [v  - vt - c3*(var'*var) - [p1, p2]*var_lift + p2u];

F = [ 
     coefficients(r - Y*N, var_aug) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     sos(v - var'*var),...
     cv == cvt,...
     ];

opts = sdpsettings('solver','mosek','verbose',0);

sol = solvesos(F, [], opts, [c1;c2;c2u;cv;cvt;coefficients(Y,var_aug)])

%% extract solution Thm1
syms z1 z2 z3 z4
var = [z1;z2];
var_aug = [z1;z2;z3;z4];

cv = value(cv);   %  cv = cv.*(cv>=1e-6 | cv<=-1e-6);    % clean small values
vv = monomials(var, 1:dv);
V = cv'*vv;
V = matlabFunction(V);

% c1  = value(c1);
% c2  = value(c2);  %  c2  = c2.*(c2>=1e-6 | c2<=-1e-6);
% c2u = value(c2u); %  c2u = c2u.*(c2u>=1e-6 | c2u<=-1e-6);
% c2u(1) = 0;
% 
% vp = monomials(var_aug, 0:dp);
% P1 = c1'*vp;
% P2 = c2'*vp;
% P2U = c2u'*vp;
% 
% P1 = matlabFunction(P1);
% P2 = matlabFunction(P2);
% P2U = matlabFunction(P2U);

VF = monomials(var,1:df);
VF = matlabFunction(VF);

% get system in syms
[f,g] = getSystem(eg,var); 
FF = matlabFunction(f);

% % check if V is sos
figure()
fc = fcontour(V,[-2 2 -2 2]);
fc.LevelList = [0:20];
% colorbar
axis square
hold on

%% 0301 Thm 2

% initial point 
x = [2;1];
% x = [0.048; 0.7737];
UU = [];
XX = x;
VV = V(x(1),x(2));

for j = 1:40
    fprintf('%d-th iteration \n',j)
    % yalmip setup
    yalmip('clear')
    sdpvar u

    % vf
    vf = VF(x(1),x(2));

    % optimize
    F = [];
    for i = 1:length(Vert)
        xt = kron(eye(2),vf')*Vert(i,:)'+g*u;
        F = [F;
             V(xt(1),xt(2)) <= V(x(1),x(2))-c3*(x'*x);
             ];
    end
    obj = norm(u);

    opts = sdpsettings('solver','mosek','verbose', 0);
    sol = optimize(F,obj,opts);

    u = double(u);
    xnext = FF(x(1),x(2))+g*u;
    dV = V(x(1),x(2)) - V(xnext(1),xnext(2));
    
    x = xnext;
    UU = [UU, u];
    XX = [XX, xnext];
    VV = [VV, V(xnext(1),xnext(2))]; 

end

plot(XX(1,:),XX(2,:))

figure()
plot(1:length(VV),VV)
legend('lyap')
