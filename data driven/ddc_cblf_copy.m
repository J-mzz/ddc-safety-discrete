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
N = [A B; -A -B];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end
e = [eps*ones(Cons.n*Cons.T,1)+xi; eps*ones(Cons.n*Cons.T,1)-xi];

% keep nontrivial constraints
[N, e] = nontrivial_constraints(N, e);

%% Unsafe region

r0 = 1;         % radius
c0 = [1; -2];	% center

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

% vp and vf
[~,~,vf] = polynomial(var,df,1);
[~,~,vg] = polynomial(var,dg,0);

% auxiliary varible Y
Y = zeros(1, size(N,1), 'like', sdpvar);	% the polynomial matrix Y
for i = 1:size(N,1)
    Y(i) = polynomial(var_aug, dv, 0);
end

% for stable and safety set P2
r = [ - kron([p1, p2], vf'), - kron([0, p2u], vg')];
d = [v  - vt  - [p1, p2]*var_lift - 1e-4*var'*var];

F = [ 
     coefficients(r - Y*N, var_aug) == 0,...
     sos(-Y*e + d),...
     sos(Y),...
     sos(v - var'*var),...
     cv == cvt,...
     ];

opts = sdpsettings('solver','mosek','verbose',0);

sol = solvesos(F, [], opts, [c1;c2;c2u;cv;cvt;coefficients(Y,var_aug)])

%% extract solution

syms z1 z2 z3 z4
var = [z1;z2];
var_aug = [z1;z2;z3;z4];

cv = value(cv);     cv = cv.*(cv>=1e-6 | cv<=-1e-6);    % clean small values
c1 = value(c1);     c1 = c1.*(c1>=1e-9 | c1<=-1e-9);
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
P1 = c1'*v1; 
P2 = c2'*v2; 
P2U = c2u'*v2u;

U = P2U/P2;     % extract u = p2u / p2

V = matlabFunction(V);
P1 = matlabFunction(P1);
P2 = matlabFunction(P2);
P2U = matlabFunction(P2U);
U = matlabFunction(U);

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

%% test trajectory with p(x, tilde x)

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

        for i = 1:150
            next = FF(curr(1),curr(2)) + g*u;

            g1 = V(curr(1),curr(2)) - V(next(1),next(2)) - tol;%*0.9^i;
            P = msdp(min(g0), [g1 >= 0], 3);

            [status,obj] = msol(P);

            if status ~= 1
        %        checkV(traj, V)
                warning('Numerical issue occurs at i=%d iteration', i)
        %        traj(:,end)=[];
                break
            end

            ulog = [ulog, double(u)];
            traj = [traj, double(next)];
            Vlog = [Vlog, V(double(next(1)),double(next(2)))];

            curr = double(next);

        end

        plot(traj(1,:),traj(2,:),'--')
        plot(traj(1,end),traj(2,end),'sr')


        theta = (1:30)/30*pi;
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

%% auxiliary functions 

% check V(traj)
function checkV(traj, V)
    for j = 1:size(traj,2)
        V(traj(1,j), traj(2,j))
    end
end