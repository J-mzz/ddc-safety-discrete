function x = getData(eg,tspan,x_init,eps)
% Starting from given initial points, generate and randomly sample the
% trajectories, then add noise. The input is a random number between [-1,1]
%
% System dynamics: x_dot = f(x) + g(x)u(x) + noise

global x_glo
global index

x = [];
x_glo = zeros();

for i = 1:size(x_init,2)
    index = 1;
    x0 = x_init(:,i);
    
    [~, ~] = ode45(@(t,x) sysDyn(t,x,eps,eg), tspan, x0);

    x = [x;x_glo];

end

clearvars -global

end



%%
function dxdt = sysDyn(t,x,eps,eg)

global x_glo
global index

% u = 2*rand(1,1) - 1;
u = -1;

x_glo(index, 1) = u(1);
x_glo(index, 2) = t;
x_glo(index, 3) = x(1);
x_glo(index, 4) = x(2);

noise = eps*(2*rand(2,1)-1);


switch eg
    case 1
        A = [-0.6409    0.8962
            0.2408   -0.6790];
        g = [0; 1];
        
        dxdt = A*[x(1);x(2)] + g*u + noise;
        
    case 2
        g = [0; 1];
        
        dxdt = [x(2); -x(1) + 1/3*x(1)^3 - x(2)] + g*u + noise;
        
    case 3
        A = [-0.6409    0.8962
            0.2408   -0.6790];
        g = [0; 1];
        
        dxdt = -A*[x(1);x(2)] + g*u + noise;

end

% dxdt = [x(2) - x(1)^3 + x(1)^2; 0] + [1;2]*u + noise;

x_glo(index, 5) = dxdt(1);
x_glo(index, 6) = dxdt(2);

index = index + 1;

end