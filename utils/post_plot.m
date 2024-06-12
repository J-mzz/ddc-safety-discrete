close all;clear;clc
%% unsafe trajectory
load unsafe_2dconvex.mat
figure()
unsafe_traj = plot(XX(1,:),XX(2,:),'--*','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
hold on

%% safe trajectory
load safe_2dconvex.mat

fc = fcontour(V,[-2 2 -2 2]);
fc.LevelList = [0:20];
% colorbar
axis square

safe_traj = plot(XX(1,:),XX(2,:),'--*','LineWidth',1.5,'Color',[0 0.4470 0.7410]);

if isLinearH
    X_unsafe = -2:.1:2;
    Y_unsafe = -1*ones(1,length(X_unsafe));
    unsafe_reg = plot(X_unsafe,Y_unsafe,'r','LineWidth',1.5);
else
    theta = (0:60)/30*pi;
    X_unsafe = c0(1) + sqrt(r0)*cos(theta);
    Y_unsafe = c0(2) + sqrt(r0)*sin(theta);
    unsafe_reg = plot(X_unsafe,Y_unsafe,'r','LineWidth',1.5);
end

xlabel('\boldmath${x_1}$','Interpreter','latex')
ylabel('\boldmath${x_2}$','Interpreter','latex')

% xlim([-2 2])
% ylim([-2 2])
xlim([-.5 2])
ylim([-1 1.5])
legend([safe_traj, unsafe_traj, unsafe_reg],'Safe traj','Unsafe traj','Unsafe reg','FontSize',14) 


% %% safe Lyapunov trace
% figure()
% subplot(2,1,1);
% safe_L = plot(0:length(VV)-1,VV,'LineWidth',1.5,'Color',[0 0.4470 0.7410]);
% hold on
% 
% 
% %% unsafe Lyapunov trace
% load unsafe_2dnonconvex.mat
% unsafe_L = plot(0:length(VV)-1,VV,'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% 
% % xlabel('\boldmath${k}$','Interpreter','latex')
% ylabel('\boldmath${V(x_k)}$','Interpreter','latex')
% 
% legend([safe_L, unsafe_L],'Safe trace','Unsafe trace','FontWeight','bold','FontSize',12)
% hold off
% 
% 
% 
% %% unsafe control input
% subplot(2,1,2);
% unsafe_con = plot(0:length(UU)-1,UU,'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% hold on
% 
% 
% %% safe control input
% load safe_2dnonconvex.mat
% 
% safe_con = plot(0:length(UU)-1,UU,'LineWidth',1.5,'Color',[0 0.4470 0.7410]);
% xlabel('\boldmath${k}$','Interpreter','latex')
% ylabel('\boldmath${u(x_k)}$','Interpreter','latex')
% 
% legend([safe_con, unsafe_con],'Safe control','Unsafe control','FontWeight','bold','FontSize',12)
% hold off
% 


%% safe Lyapunov trace
figure()
yyaxis left
safe_L = plot(0:length(VV)-1,VV,'-*','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
hold on
axis square

%% unsafe Lyapunov trace
load unsafe_2dconvex.mat
unsafe_L = plot(0:length(VV)-1,VV,'-*','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);

% xlabel('\boldmath${k}$','Interpreter','latex')
ylabel('\boldmath${V(x_k)}$','Interpreter','latex')

% legend([safe_L, unsafe_L],'Safe trace','Unsafe trace','FontWeight','bold','FontSize',12)



%% unsafe control input
yyaxis right
unsafe_con = plot(0:length(UU)-1,UU,':.','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);


%% safe control input
load safe_2dconvex.mat

safe_con = plot(0:length(UU)-1,UU,':.','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
xlabel('\boldmath${k}$','Interpreter','latex')
ylabel('\boldmath${u(x_k)}$','Interpreter','latex')

legend([safe_L, unsafe_L, safe_con, unsafe_con],'Safe V','Unsafe V','Safe u','Unsafe u','FontWeight','bold','FontSize',12)
hold off

