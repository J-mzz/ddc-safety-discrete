close all;clear;clc
%% unsafe trajectory
load unsafe_3dconvex.mat
figure()
unsafe_traj = plot3(XX(1,:),XX(2,:),XX(3,:),'--*','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
hold on
grid on
%% safe trajectory
load safe_3dconvex.mat

safe_traj = plot3(XX(1,:),XX(2,:),XX(3,:),'--*','LineWidth',1.5,'Color',[0 0.4470 0.7410]);

unsafe1 = [1  1  1  1];
unsafe2 = [-2 2  2 -2];
unsafe3 = [-2 -2 2  2];
unsafe_reg = fill3(unsafe1,unsafe2,unsafe3,'r');

xlabel('\boldmath${x_1}$','Interpreter','latex')
ylabel('\boldmath${x_2}$','Interpreter','latex')
zlabel('\boldmath${x_3}$','Interpreter','latex')

xlim([-1;2])
caz =-15;
cel = 54;
view(caz,cel)
hold off

legend([safe_traj, unsafe_traj, unsafe_reg],'Safe traj','Unsafe traj','Unsafe reg','FontSize',14) 


% %% safe Lyapunov trace
% figure()
% subplot(2,1,1);
% safe_L = plot(0:length(VV)-1,VV,'LineWidth',1.5,'Color',[0 0.4470 0.7410]);
% hold on
% 
% 
% %% unsafe Lyapunov trace
% load unsafe_3dconvex.mat
% unsafe_L = plot(0:length(VV)-1,VV,'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
% 
% xlabel('\boldmath${k}$','Interpreter','latex')
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
% load safe_3dconvex.mat
% 
% safe_con = plot(0:length(UU)-1,UU,'LineWidth',1.5,'Color',[0 0.4470 0.7410]);
% xlabel('\boldmath${k}$','Interpreter','latex')
% ylabel('\boldmath${u(x_k)}$','Interpreter','latex')
% 
% legend([safe_con, unsafe_con],'Safe control','Unsafe control','FontWeight','bold','FontSize',12)
% hold off


%% safe Lyapunov trace
figure()
yyaxis left
safe_L = plot(0:length(VV)-1,VV,'-*','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
hold on
axis square

%% unsafe Lyapunov trace
load unsafe_3dconvex.mat
unsafe_L = plot(0:length(VV)-1,VV,'-*','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);

% xlabel('\boldmath${k}$','Interpreter','latex')
ylabel('\boldmath${V(x_k)}$','Interpreter','latex')

% legend([safe_L, unsafe_L],'Safe trace','Unsafe trace','FontWeight','bold','FontSize',12)



%% unsafe control input
yyaxis right
unsafe_con = plot(0:length(UU)-1,UU,':.','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);


%% safe control input
load safe_3dconvex.mat

safe_con = plot(0:length(UU)-1,UU,':.','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
xlabel('\boldmath${k}$','Interpreter','latex')
ylabel('\boldmath${u(x_k)}$','Interpreter','latex')

legend([safe_L, unsafe_L, safe_con, unsafe_con],'Safe V','Unsafe V','Safe u','Unsafe u','FontWeight','bold','FontSize',12)
hold off


