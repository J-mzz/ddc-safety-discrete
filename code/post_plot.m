close all;clear;clc

load data\linear_safety.mat

figure()
fc = fcontour(V,[-2 2 -2 2]);
fc.LevelList = [0:20];
% colorbar
axis square
hold on

plot(XX(1,:),XX(2,:),'LineWidth',1.5)

if isLinearH
    X_unsafe = -2:.1:2;
    Y_unsafe = -1*ones(1,length(X_unsafe));
    plot(X_unsafe,Y_unsafe,'r','LineWidth',1.5)
else
    theta = (0:60)/30*pi;
    X_unsafe = c0(1) + sqrt(r0)*cos(theta);
    Y_unsafe = c0(2) + sqrt(r0)*sin(theta);
    plot(X_unsafe,Y_unsafe,'r','LineWidth',1.5)
    ylim([-2 2])
    xlim([-2 2])
    axis square
end

xlabel('\boldmath${x_1}$','Interpreter','latex')
ylabel('\boldmath${x_2}$','Interpreter','latex')
hold off

figure()
plot(1:length(VV),VV,'LineWidth',1.5)
xlabel('\boldmath${k}$','Interpreter','latex')
ylabel('\boldmath${V(x_k)}$','Interpreter','latex')
legend('lyap','FontWeight','bold','FontSize',12)