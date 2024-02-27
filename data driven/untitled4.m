yalmip('clear')

mpol x3 x4
var = [x3;x4];              % \tilde x

% multiplier p
p1  = 0;%P1(x3,x4); % if error, let p1 = 0
p2  = P2(x3,x4);
p2u = P2U(x3,x4);

% lyapunov function v
vt = V(x3,x4);

% vp and vf
vf = VF(x(1),x(2));
vg = VG;

% new
Y = {};
for i = 1:size(N,1)
    Y{i} = YY{i}(x3,x4);
end

% for stable set P2
r = - [p1*vf', p2*vf'];
d = alpha  - V(x3,x4) - [p1, p2]*var + p2u - S(x3,x4);

Ye = 0;
for i = 1:length(Y)
    Ye = Ye+Y{i}*e(i);
end

F = [ 
     -Ye + d >= 0,...
     ];

obj = ((x3-1)*(x3-1)+(x4-2)*(x4-2));
P = msdp(min(obj),F,4)

mset('yalmip',true)
mset(sdpsettings('solver','mosek'))

[status,obj] = msol(P)