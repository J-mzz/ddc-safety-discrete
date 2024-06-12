close all;clear;clc
yalmip('clear')

rng(1)

%% system and data generation system with dim 1

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
[A,xi] = getCons2_new(X,Cons);    % with stronger a priori

% check enough data
N = [A; -A];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];


% find vertices
[N, e] = nontrivial_constraints(N, e);
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

%% Unsafe region
isLinearH = 0;
isConvex = 0;

if isLinearH
    h  = x2+1;
    ht = x4+1;
    
else
    r0 = 1;%0.25;      % radius
    c0 = [1;-1];%[1; -1.5];	% center
    if isConvex
        r0 = 4;         % radius
        c0 = [0; 1];	% center
    end
    
    h  = -r0 + (x1-c0(1))^2 + (x2-c0(2))^2; % unsafe: h(x)  < 0
    ht = -r0 + (x3-c0(1))^2 + (x4-c0(2))^2; % unsafe: h(xt) < 0
    if isConvex
        h = -h;
        ht = -ht;
    end
end

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

% multipler s: sos( H(xt)-s1*H(x) ) 
[s,cs] = polynomial(var_aug, dp, 0);

% vp and vf
[~,~,vf] = polynomial(var,df,1);
vf = vf([1 2 3]);

% auxiliary varible Y
Y = zeros(2, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(1,i) = polynomial(var_aug, dv, 0);
    Y(2,i) = polynomial(var_aug, dv, 0);
end

% for stable and safety set P2
c3 = 1e-2;
r = [- kron([p1, p2], vf');
     - kron([p1, p2], vf')]; 
d = [v  - vt - c3*(var'*var) - [p1, p2]*var_lift + p2u;
     ht - s*h                - [p1, p2]*var_lift + p2u];

F = [ 
     coefficients(r - Y*N, var_aug) == 0,...
     (sos(-Y*e + d)):'tagged',...
     sos(Y),...
     sos(v - var'*var),...
     cv == cvt,...
     ];

opts = sdpsettings('solver','mosek','verbose',0);

[sol,sol_u,sol_Q] = solvesos(F, [], opts, [c1;c2;c2u;cv;cvt;cs;coefficients(Y,var_aug)]);
sol

%% extract solution Thm1
syms z1 z2 z3 z4
var = [z1;z2];
var_aug = [z1;z2;z3;z4];

cv = value(cv);   %  cv = cv.*(cv>=1e-6 | cv<=-1e-6);    % clean small values
vv = monomials(var, 1:dv);
V = cv'*vv;
% check V sos
% findsos(V)
V = matlabFunction(V);

if isLinearH
    H = z2+1;
else
    H = -r0 + (z1-c0(1))^2 + (z2-c0(2))^2; % unsafe: h(x)  < 0
    if isConvex
        H = -H;
    end
end
H = matlabFunction(H);

VF = monomials(var,1:df);
VF = VF([1 2 3]);
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

%% 0301 Thm 2, use gloptipoly due to the safety constraint

% initial point 
x = [2;1];
UU = [];
XX = x;
VV = V(x(1),x(2));

for j = 1:10
    fprintf('%d-th iteration \n',j)
    % yalmip setup
    yalmip('clear')
    
    vf = VF(x(1),x(2));
    
    %% ??? how to solve for u
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Attemp 1
    ub = 2;
    lb = -2;
    
    deg = 4;
    
    u = sdpvar(1);
    f = sdpvar(6,1);
    
    [p,cp] = polynomial(u, deg, 0);
    
    [s01,c01] = polynomial(f, deg-1, 0);
    [s02,c02] = polynomial(f, deg-1, 0);
    
    [s11,c11] = polynomial(f, deg-1, 0);
    [s12,c12] = polynomial(f, deg-1, 0);

    [s21,c21] = polynomial(f, deg-1, 0);
    [s22,c22] = polynomial(f, deg-1, 0);

    [s31,c31] = polynomial(f, deg-1, 0);
    [s32,c32] = polynomial(f, deg-1, 0);

    [s41,c41] = polynomial(f, deg-1, 0);
    [s42,c42] = polynomial(f, deg-1, 0);

    [s51,c51] = polynomial(f, deg-1, 0);
    [s52,c52] = polynomial(f, deg-1, 0);

    [s61,c61] = polynomial(f, deg-1, 0);
    [s62,c62] = polynomial(f, deg-1, 0);
    
    % objective
    obj = 0;
    for i = 1:length(cp)
        integral = 1/i*(ub^i - lb^i);
        obj = obj + cp(i)*integral;
    end
    obj = -obj;
    
    xk1 = vf'*f(1:3) + g(1)*u;
    xk2 = vf'*f(4:6) + g(2)*u;
    
    Delta = V(x(1),x(2)) - V(xk1,xk2) - c3*(x'*x);
    
    h1 = Delta - p;
    h2 = H(xk1,xk2) - p;
    
    % constraints
    F = [   sos(s01); sos(s02); ...
            sos(s11); sos(s12); ...
            sos(s21); sos(s22); ...
            sos(s31); sos(s32); ...
            sos(s41); sos(s42); ...
            sos(s51); sos(s52); ...
            sos(s61); sos(s62); ];
        
    F = [F; sos(h1 - s01*(u   - lb)     + s02*(u   - ub)... 
                   - s11*(f(1)-vmin(1)) + s12*(f(1)-vmax(1))... 
                   - s21*(f(2)-vmin(2)) + s22*(f(2)-vmax(2))... 
                   - s31*(f(3)-vmin(3)) + s32*(f(3)-vmax(3))... 
                   - s41*(f(4)-vmin(4)) + s42*(f(4)-vmax(4))... 
                   - s51*(f(5)-vmin(5)) + s52*(f(5)-vmax(5))... 
                   - s61*(f(6)-vmin(6)) + s62*(f(6)-vmax(6)) )];
    F = [F; sos(h2 - s01*(u   - lb)     + s02*(u   - ub)... 
                   - s11*(f(1)-vmin(1)) + s12*(f(1)-vmax(1))... 
                   - s21*(f(2)-vmin(2)) + s22*(f(2)-vmax(2))... 
                   - s31*(f(3)-vmin(3)) + s32*(f(3)-vmax(3))... 
                   - s41*(f(4)-vmin(4)) + s42*(f(4)-vmax(4))... 
                   - s51*(f(5)-vmin(5)) + s52*(f(5)-vmax(5))... 
                   - s61*(f(6)-vmin(6)) + s62*(f(6)-vmax(6)) )];

    opts = sdpsettings('solver','mosek','verbose',0);
    coeff = [cp;c01;c02;c11;c12;c21;c22;c31;c32;c41;c42;c51;c52;c61;c62];
    sol = solvesos(F, obj, opts, coeff)
%     if sol.problem
%         error('error')
%     end

    
    % extract p
    syms u
    cp = value(cp);
    vp = monomials(u, 0:deg);
    P = cp'*vp;
    P = matlabFunction(P);
    
    t = lb:.01:ub;
    for i = 1:length(t)
        Pu(i) = P(t(i));
    end
    
    % find u with 
    [~,ind] = max(Pu);
    u = t(ind);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Attemp 1
    
    % update system dynamics
    xnext = FF(x(1),x(2))+g*u;
    dV = V(x(1),x(2)) - V(xnext(1),xnext(2));
    
    x = xnext;
    UU = [UU, u];
    XX = [XX, xnext];
    VV = [VV, V(xnext(1),xnext(2))]; 

end

plot(XX(1,:),XX(2,:),'--*')

if isLinearH
    X_unsafe = -2:.1:2;
    Y_unsafe = -1*ones(1,length(X_unsafe));
    plot(X_unsafe,Y_unsafe,'r')
else
    theta = (0:60)/30*pi;
    X_unsafe = c0(1) + sqrt(r0)*cos(theta);
    Y_unsafe = c0(2) + sqrt(r0)*sin(theta);
    plot(X_unsafe,Y_unsafe,'r')
    ylim([-2 2])
    xlim([-2 2])
    axis square
end
xlabel('x_1')
ylabel('x_2')
hold off

% figure()
% plot(1:length(VV),VV)
% xlabel('k')
% ylabel('V(x_k)')
% legend('lyap')
