close all;clear;clc
yalmip('clear')

rng(1)

%% system and data generation

% yalmip setup
sdpvar x1 x2 x3 x4
var = [x1;x2];              % x
var_lift = [x3;x4];         % \tilde{x}: x3 = next(x1), x4 = next(x2)
var_aug = [var; var_lift];  % [x; \tilde{x}]

% define DT system: x+ = f(x) + g(x)*u(x) + noise
eg = 1;                     % 1 (linear) or 2 (nonlinear)
[f,g] = getSystem(eg,var); 

% setup
Cons = setCons(eg);

% noise level
eps = 1e-1;                    % noise level for robust data driven
epsw = 0;                   % noise level of disturbance

% generate data, with input u \in [-1,1]
X = getDataRND(eg,eps,Cons.T);

% consistency set
[A,B,xi] = getCons(X,Cons);

% check enough data
N = [A; -A];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];

% keep nontrivial constraints
[N, e] = nontrivial_constraints(N, e);

%% Unsafe region

r0 = 1;         % radius
c0 = [1; -2];	% center

h  = -r0 + (x1-c0(1))^2 + (x2-c0(2))^2; % unsafe: h(x)  < 0
ht = -r0 + (x3-c0(1))^2 + (x4-c0(2))^2; % unsafe: h(xt) < 0


% h  = -1 + x1*x2*(-9.811350740483604e-2)+x1^2*3.663287817314286+x2^2*2.613869172425404;
% ht = -1 + x3*x4*(-9.811350740483604e-2)+x3^2*3.663287817314286+x4^2*2.613869172425404;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try to solve
n   = Cons.n;
T   = Cons.T;
df  = Cons.df;
dg  = Cons.dg;

%% change the value of FLAG to switch from p(x) and p(x,\tilde{x})

FLAG = 1;   % FLAG = 0 -->  p(x)
            % FLAG = 1 -->  p(x,\tilde{x})

%% SOS set up

dp = 2;     % degree of vecter polynomial p = [p1; p2]
dv = 4;     % degree of lyapunov function v

% multiplier p
if ~FLAG
    % p(x)
    [p1,c1] =   polynomial(var, dp, 0);
    [p2,c2] =   polynomial(var, dp, 0);
    [p2u,c2u] = polynomial(var, dp, 0);     % p2u = p2*u
else
    % p(x, \tilde{x})
    [p1,c1] =   polynomial(var_aug, dp, 0);
    [p2,c2] =   polynomial(var_aug, dp, 0);
    [p2u,c2u] = polynomial(var_aug, dp, 0); % p2u = p2*u
end

% lyapunov function v
[v, cv] =  polynomial(var, dv, 1);          % 1 -> without constant
[vt,cvt] = polynomial(var_lift, dv, 1);

% multipler s: sos( H(xt)-s1*H(x) ) 
[s,cs] = polynomial(var_aug, dp, 0);

% vp and vf
[~,~,vf] = polynomial(var,df,1);
[~,~,vg] = polynomial(var,dg,0);

% new
[~,~,vY] = polynomial(var_aug, dv, 0);

% auxiliary varible Y
Y = zeros(2, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(1,i) = polynomial(var_aug, dv, 0);
    Y(2,i) = polynomial(var_aug, dv, 0);
end

% for stable and safety set P2
r = [ - kron([p1, p2], vf');
      - kron([p1, p2], vf') ];
d1 = v  - vt  - [p1, p2]*var_lift + p2u -1e-4 * var'*var; 
d2 = ht - s*h - [p1, p2]*var_lift + p2u;
d = [d1;d2];

F = [ 
     coefficients(r - Y*N, var_aug) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     sos(s),...
     sos(v - var'*var),...
     cv == cvt,...
     ];

opts = sdpsettings('solver','mosek','verbose',0);

sol = solvesos(F, [], opts, [c1;c2;c2u;cv;cvt;cs;coefficients(Y,var_aug)])

%% extract solution

cY = value(coefficients(Y,var_aug));
cY = cY.*(cY>=1e-6 | cY<=-1e-6);

syms z1 z2 z3 z4
var = [z1;z2];
var_aug = [z1;z2;z3;z4];

cv = value(cv);     cv = cv.*(cv>=1e-6 | cv<=-1e-6);    % clean small values
% c1 = value(c1);     c1 = c1.*(c1>=1e-9 | c1<=-1e-9);
c2 = value(c2);     c2 = c2.*(c2>=1e-9 | c2<=-1e-9);
c2u = value(c2u);   c2u = c2u.*(c2u>=1e-6 | c2u<=-1e-6);

vv = monomials(var, 1:dv);
if ~FLAG
    % p(x)
    v1 =  monomials(var, 0:dp);
    v2 =  monomials(var, 0:dp);
    v2u = monomials(var, 0:dp);
else
    % p(x, \tilde{x})
    v1 =  monomials(var_aug, 0:dp); 
    v2 =  monomials(var_aug, 0:dp);
    v2u = monomials(var_aug, 0:dp);
end

V = cv'*vv;
% P1 = c1'*v1; 
P2 = c2'*v2; 
P2U = c2u'*v2u;

U = P2U/P2;     % extract u = p2u / p2

vY = monomials(var_aug, 0:dv);
nY = length(vY);
for i = 1:size(N,1)
    YY{1,i} = cY((1:nY)+nY*(2*i-2))'*vY;
    YY{2,i} = cY((1:nY)+nY*(2*i-1))'*vY;
    
    YY{1,i} = matlabFunction(YY{1,i});
    YY{2,i} = matlabFunction(YY{2,i});
end

V = matlabFunction(V);
% P1 = matlabFunction(P1);
P2 = matlabFunction(P2);
P2U = matlabFunction(P2U);
U = matlabFunction(U);

H = -r0 + (z1-c0(1))^2 + (z2-c0(2))^2; % unsafe: h(x)  < 0
H = matlabFunction(H);

VF = monomials(var,1:df);
VG = monomials(var,0:dg);

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


%% Online method 4. Error: Missing non-negativity bounds in GP formulation
%{
% online optimization at point x
x = [2;1];

% yalmip setup
yalmip('clear')
sdpvar s1 s2
sdpvar a
sdpvar x3 x4
var = [x3;x4];              % \tilde x

% multiplier p
% p1  = P1(x(1),x(2),x3,x4);
p1 = 0;
p2  = P2(x(1),x(2),x3,x4);
p2u = P2U(x(1),x(2),x3,x4);

% lyapunov function v
vt = V(x3,x4);

% lyapunov function v
ht = H(x3,x4);

% vp and vf
vf = VF(x(1),x(2));
vg = VG;

% new
Y = zeros(2, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(1,i) = YY{1,i}(x(1),x(2),x3,x4);
    Y(2,i) = YY{2,i}(x(1),x(2),x3,x4);
end

% for stable and safety set P2
r = [ - kron([p1, p2], vf');
      - kron([p1, p2], vf') ];
d1 = a  - vt - [p1, p2]*var + p2u - s1;
d2 = ht      - [p1, p2]*var + p2u - s2;

F = [ 
%      r == Y*N,...
%      -Y(1,:)*e + d1 >= 0,...
     -Y(2,:)*e + d2 >= 0,...
     s1 >= 0,...
     s2 >= 0,...
     a >= V(x(1),x(2))
%      Y >= 0,...
     ];

opts = sdpsettings('solver','mosek','verbose',1);

sol = optimize(F, a, [])

u = P2U(x(1),x(2),value(x3),value(x4))/P2(x(1),x(2),value(x3),value(x4));


%}

%% Online method 3
%
x = [2;1];
%
yalmip('clear')

% yalmip setup
sdpvar a
sdpvar x3 x4
var = [x3;x4];              % tilde x

% multiplier p
[p1,c1] =   polynomial(var, dp, 0);
[p2,c2] =   polynomial(var, dp, 0);
[p2u,c2u] = polynomial(var, dp, 0);     % p2u = p2*u

% multipler s
[s1,cs1] = polynomial(var, dp, 0);
[s2,cs2] = polynomial(var, dp, 0);

% lyapunov v, barrier h
[v, cv] =  polynomial(var, dv, 1);          % 1 -> without constant
ht  = -r0 + (x3-c0(1))^2 + (x4-c0(2))^2; % unsafe: h(x)  < 0

% vp and vf
vf = VF(x(1),x(2));
vg = VG;

% auxiliary varible Y
Y = zeros(2, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(1,i) = polynomial(var, dv, 0);
    Y(2,i) = polynomial(var, dv, 0);
end

% for stable and safety set P2
% r = [ - kron([p1, p2], vf'), - kron([0, p2u], vg');
%       - kron([p1, p2], vf'), - kron([0, p2u], vg') ];
% d = [a        - V(x3,x4) - [p1, p2]*var - s1;
%      H(x3,x4) - s2       - [p1, p2]*var  ];

r = [ - kron([p1, p2], vf');
      - kron([p1, p2], vf')];
d1 = a  - V(x3,x4)   - [p1, p2]*var +p2u - s2;
d2 = ht - s1         - [p1, p2]*var +p2u;
d  = [d1;d2];

F = [ 
     coefficients(r - Y*N, var) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     a >= V(x(1),x(2)),...
     ];

opts = sdpsettings('solver','mosek','verbose',1);

sol = solvesos(F, a, [], [c1;c2;c2u;cv;cs1;cs2;coefficients(Y,var);a])
%}
%% Online method 1
%{
yalmip('clear')
mset('yalmip',true)
mset(sdpsettings('solver','mosek','verbose',1))
N1 = size(N,1);

mpol alpha
mpol x1 x2
mpol p1 p2 p2u
mpol y1 1 104
mpol y2 1 104
mpol s

xk = [2;1];

% vp and vf
vf = monomials([z1;z2],0:df);
VF = matlabFunction(vf);

% for stable and safety set P2
r  = [-p1*VF(xk(1),xk(2))', -p2*VF(xk(1),xk(2))'];
d1 = alpha    - V(xk(1),xk(2))   - [p1, p2]*[x1;x2] + p2u;
d2 = H(x1,x2) - s*H(xk(1),xk(2)) - [p1, p2]*[x1;x2] + p2u;

F = [ 
     r == y1*N,...
     r == y2*N,...
     d1 >= y1*e,...
     d2 >= y2*e,...
     y1 >= 0,...
     y2 >= 0,...
     s >= 0,...
     ];

P = msdp(min(alpha), F, 2);
[status,obj] = msol(P);
%}
%% Online method 2
%{
yalmip('clear')

sdpvar alpha
sdpvar x1 x2
sdpvar p1 p2 p2u
y1 = sdpvar(1, size(N,1));
y2 = sdpvar(1, size(N,1));
sdpvar s
sdpvar p1x p2x

xk = [2;1];

% vp and vf
vf = monomials([z1;z2],0:df);
VF = matlabFunction(vf);

% for stable and safety set P2
% r  = -kron([p1,p2], VF(xk(1),xk(2))');
r = -[p1*VF(xk(1),xk(2))', p2*VF(xk(1),xk(2))'];
d1 = alpha    - V(xk(1),xk(2))   - p1x - p2x + p2u;
d2 = H(x1,x2) - s*H(xk(1),xk(2)) - p1x - p2x + p2u;

F = [ 
     r - y1*N == 0,...
     r - y2*N == 0,...
     d1 >= y1*e,...
     d2 >= y2*e,...
     y1 >= 0,...
     y2 >= 0,...
     s >= 0,...
     ];

opts = sdpsettings('solver','mosek','verbose',1);

sol = optimize(F, alpha, opts)
%}
%% WRONG WAY!!! test trajectory with p(x, tilde x)
%{
switch eg

%%  test trajectory using yalmip
% NOTE: this only works when V is quadratic of x!
% if not, use the next block

    case 1
        mpol u
        tol = 1e-1;
        if eps == 0
            tol = 0;
        end

        curr = [2; 1];
        traj = curr;
        ulog = [];
        Vlog = V(curr(1),curr(2));
        plot(curr(1),curr(2),'or')

        % objective
        g0 = u^2;

        % setting
        mset('yalmip',true)
        mset(sdpsettings('solver','mosek','verbose',0))
        flag = 0;
        for i = 1:150
            next = FF(curr(1),curr(2)) + g*u;

            g1 = V(curr(1),curr(2)) - V(next(1),next(2)) - tol;%*0.9^i;
            g2 = H(next(1),next(2));
            P = msdp(min(g0), [g1 >= 0; g2 >= 0], 3);

            [status,obj] = msol(P);

            if status ~= 1
%                 checkV(traj, V)
                warning('Numerical issue occurs at i=%d iteration', i)
%                 traj(:,end)=[];
                flag = 1;
%                 break
            end
            if flag
                P = msdp(min(g0), [g1 >= 0], 3);
                [status,obj] = msol(P);
                if status ~= 1
%                 checkV(traj, V)
                    warning('Numerical issue occurs at i=%d iteration', i)
%                 traj(:,end)=[];
                break
                end
            end
            ulog = [ulog, double(u)];
            traj = [traj, double(next)];
            Vlog = [Vlog, V(double(next(1)),double(next(2)))];

            curr = double(next);

        end

        plot(traj(1,:),traj(2,:),'--')
        plot(traj(1,end),traj(2,end),'sr')


        theta = (0:30)/30*pi;
        X_unsafe = c0(1) + sqrt(r0)*cos(theta);
        Y_unsafe = c0(2) + sqrt(r0)*sin(theta);
        plot(X_unsafe,Y_unsafe)
        hold off

        figure()
        plot(1:length(Vlog),Vlog)



%% test trajectory using gloptipoly
% use gloptipoly if yalmip does not work... 
    case 2 %%%% DO NOT USE RIGHT NOW!!!
        mpol u
        tol = 1e-1;
        if eps==0
            tol = 0;
        end

        curr = [1; -2]; %[2; 1];
        traj = curr;
        ulog = [];
        Vlog = V(curr(1),curr(2));
        plot(curr(1),curr(2),'or')
        
        % objective
        g0 = u^2;
        
        % setting
        mset('yalmip',true)
        mset(sdpsettings('solver','mosek','verbose',0))
        
        for i = 1:150
            next = FF(curr(1),curr(2)) + g*u;
        
            g1 = V(curr(1),curr(2)) - V(next(1),next(2)) - tol*0.9^i;
            P = msdp(min(g0), g1 >= 0, 4);
        
            [status,obj] = msol(P);
            
            if status ~= 1
%                 checkV(traj, V)
                warning('Numerical issue occurs at i=%d iteration', i)
%                 traj(:,end)=[];
                break
            end
        
            ulog = [ulog, double(u)];
            traj = [traj, double(next)];
            Vlog = [Vlog, V(double(next(1)),double(next(2)))];
            
            curr = double(next);
        
        end
        
        plot(traj(1,:),traj(2,:),'--')
        plot(traj(1,end),traj(2,end),'sr')
        hold off
        
        figure()
        plot(1:length(Vlog),Vlog)
        
        
        % % test traj with no control
        % curr = [1; -2];
        % traj = curr;
        % for i = 1:10
        %     next = FF(curr(1),curr(2));
        %     traj = [traj, next];
        %     curr = next;
        % end
        % plot(traj(1,:),traj(2,:),'--')
end
%}
%% auxiliary functions 

% check V(traj)
function checkV(traj, V)
    for j = 1:size(traj,2)
        V(traj(1,j), traj(2,j))
    end
end