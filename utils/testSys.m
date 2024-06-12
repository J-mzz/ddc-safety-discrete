 function DATA = testSys(Region,f,g,u,epsw)
%% plot closed loop phase protrait

r0 = Region.r0;
c0 = Region.c0;     % initial

n = 10; % # of samples
sample = (1:n)/n * 2*pi;
x = c0(1) + sqrt(r0) * sin(sample);
y = c0(2) + sqrt(r0) * cos(sample);


FF = matlabFunction(f+g*u);
F = matlabFunction(f);
U = matlabFunction(u);
% G = matlabFunction(g);
%-------------------------Trajectory: ode45

% tspan = [0:.2:2]; % time span
% curr_ode_options =  odeset('RelTol', 1e-7, ...
%     'AbsTol', 1e-8, 'MaxStep', 1e-6);
% 
% Noise = @(t, x) epsw*(2*rand(2,1)-ones(2,1));
% FFF = @(t,x) FF(x(1),x(2)) + Noise(t,[x(1),x(2)]);
% 
% for i = 1:n
%     x0 = [0;-2.5]; %[x(9);y(9)];0.5;-3
%     X = ode15s(FFF, tspan, x0, curr_ode_options);   % u \in [-1,1]
%     
%     DataX{i} = X.y(1,:);
%     DataY{i} = X.y(2,:);
%     i
% % 	plot(X.y(1,:), X.y(2,:),'c','LineWidth',2);
%     plot(X.y(1,:), X.y(2,:),'b');
% end
% DATA.X = DataX;
% DATA.Y = DataY;
%% plot 3d trajectories
x0 = [0;-2.5];

%%
%
T = 2;
mu = 0.1;
ode_opt = odeset('RelTol', 1e-6,'AbsTol', 1e-8);

% DATA.Traj{1} = [];
% DATA.Input{1} = [];
% DATA.Time{1} = [];
% DATA.Switch{1} = [];
n = 30;

theta = (1:n)/n * 2*pi;
xxx = sqrt(r0)*sin(theta)+c0(1);
yyy = sqrt(r0)*cos(theta)+c0(2);

for i = 1:n
%     yprev = [xxx(12);yyy(12)];
    yprev = [xxx(i);yyy(i)];
%     yprev = x0;
    t_all = 0;
    switch_times = 0;

    Ydata = yprev;
    Udata = U(yprev(1),yprev(2));
    Tlog = 0;
    
    fprintf('current traj: %d \n', i)
    
    while t_all < T
        tmax_curr = exprnd(mu);
        tmax_curr = min(t_all + tmax_curr, T) - t_all;

        noise_curr = epsw*(2*rand(2,1)-1);
        
        % open loop
%         if i == 1
%             f_curr = @(t, y) F(y(1),y(2));
%         else 
%             f_curr = @(t, y) F(y(1),y(2)) + noise_curr;
%         end
        
        % close loop
%         if i == 1
%             f_curr = @(t, y) FF(y(1),y(2));
%         else 
            f_curr = @(t, y) FF(y(1),y(2)) + noise_curr;
%         end2

        [tcurr, ycurr] = ode15s(f_curr, [0, tmax_curr], yprev, ode_opt);

        Ydata = [Ydata, ycurr'];
        Udata = [Udata, U(ycurr(1),ycurr(2))];
        Tlog = [Tlog; t_all + tcurr];
        switch_times = [switch_times; t_all + tmax_curr];

        yprev = ycurr(end, :)';
        t_all = t_all + tmax_curr;
    end
    
%     if i == 1
%         plot(Ydata(1,:),Ydata(2,:),'r','LineWidth',1)
%     else
        plot(Ydata(1,2:end),Ydata(2,2:end),'b')
%     end
    
    DATA.Traj{i} = Ydata;
    DATA.Input{i} = Udata;
    DATA.Time{i} = Tlog;
    DATA.Switch{i} = switch_times;
end
%}
%-------------------------Trajectory: save with data
%{
global DATALOG
global data_index
global trajectory_index
tspan = [0:.2:2]; % time span
curr_ode_options =  odeset('RelTol', 1e-7, ...
    'AbsTol', 1e-8, 'MaxStep', 1e-6, 'Events',@myEvent);
trajectory_index = 1;
for i = 1:n
    data_index = 1;
    x0 = [0;-2.5]; %[x(9);y(9)]
    X = ode15s(@(t,x) mypoly(t,x,epsw,U(x(1),x(2))), tspan, x0, curr_ode_options);   % u \in [-1,1]
    trajectory_index = trajectory_index +1;
    data_index = 1;
%     DataX{i} = X.y(1,:);
%     DataY{i} = X.y(2,:);
    i
% 	plot(X.y(1,:), X.y(2,:),'c','LineWidth',2);
    plot(X.y(1,:), X.y(2,:),'b');
end
DATA.X = DataX;
DATA.Y = DataY;
%}
%-------------------------Trajectory: streamline
%{
% q = 5;     % bound
% k = 0.1;   % step size
% [XX,YY] = meshgrid(-q:k:q, -q:k:q);
% 
% for i = 1:size(XX,1)
%     for j = 1:size(XX,2)
% %         noise_value  = epsw*(2*rand(2,1)-ones(2,1));
% %         noise_weight = norm(epsw*(2*rand(2,1)-ones(2,1)),inf)/...
% %             U(XX(i,j),YY(i,j));
%         temp = FF(XX(i,j),YY(i,j)) + epsw*(2*rand(2,1)-ones(2,1));
%         UU(i,j) = temp(1);
%         VV(i,j) = temp(2);
%     end
% end
% % V = V + subs(u,{x1 x2}, {X Y});
% 
% % quiver(XX,YY,UU,VV, 'autoscalefactor', 10)
% h = streamline(XX,YY,UU,VV,x,y,k);
% 
% set(h,'Color','c')
% 
% axis(1*[-q q -q q])
% xlabel('x_1', 'fontsize', 13)
% ylabel('x_2', 'fontsize', 13)
% set(gca,'fontsize', 13)
%}
%-------------------------

xlim([-5 5])
ylim([-5 5])
axis square




%% plot phase portrait

plotpp(@(t,x)  F(x(1),x(2)) + g * U(x(1),x(2))) 

end


function dxdt = mypoly(t,x,epsw,u)
global DATALOG
global data_index
global trajectory_index
    u = max(min(u,1),-1);         % option 1: bound u within [-1,1]
%     if (u <= -1) || (u >= 1)      % option 2: large u is caused with
%                                               numerical issues, set to 0  
%         u = 0;
%     end
    DATALOG{trajectory_index}(data_index, 1) = u;
    DATALOG{trajectory_index}(data_index, 2) = t;
    DATALOG{trajectory_index}(data_index, 3) = x(1);
    DATALOG{trajectory_index}(data_index, 4) = x(2);
    noise = epsw*(2*rand(2,1)-1);
    dxdt = [x(2) + noise(1); -x(1)+1/3*x(1)^3-x(2) + u + noise(2)];
    DATALOG{trajectory_index}(data_index, 5) = dxdt(1);
    DATALOG{trajectory_index}(data_index, 6) = dxdt(2);
    DATALOG{trajectory_index}(data_index, 7) = noise(1);
    DATALOG{trajectory_index}(data_index, 8) = noise(2);
    data_index = data_index + 1;
end


function [value, isterminal, direction] = myEvent(t, x)
% terminate the integrator when trajectories are out of interest area
    value      = (x(1) <= -5) | (x(1) >= 5) | (x(2) <= -5) | (x(2) >= 5);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end