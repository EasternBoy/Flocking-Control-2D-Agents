% clc
% clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T  = 80;
dT = 0.5e-2;
Ns = T/dT;
N  = 19;          % Number of agents

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Agent dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nv = 4*ones(1,N);  % 4-vertices agents
po = zeros(2,4,N);

for i = 1:N
    po(:,:,i)  = [0.4 0.2; 0.4 -0.2; -0.4 -0.2; -0.4 0.2]';
end

S    = [0 -1;1 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rendezous tranjactory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rod = 5;
Fre = 1/(25*pi);
Ome = Fre*2*pi;
qod   = zeros(3,Ns);
dqod  = zeros(3,Ns);
ddqod = zeros(3,Ns);
for k = 1:1:Ns
    %% Rendezous tranjactory
    t = k*dT;
    qod(:,k)   = [Rod*cos(Ome*t+pi) + Rod;   Rod*sin(Ome*t);       -Ome*t + pi/2];
    dqod(:,k)  = [Rod*Ome*sin(Ome*t);        Rod*Ome*cos(Ome*t);   -Ome];
    ddqod(:,k) = [Rod*Ome^2*cos(Ome*t);     -Rod*Ome^2*sin(Ome*t); 0];
end

% load('q_traj.mat');
% load('z_traj.mat');
% load('up.mat');
% load('dist.mat');
% load('uphi.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf(figure(4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1);
plot_agents(q_traj,qod,po,Nv,dT,dT);
axis([-5.5   12.8   -7.5   8]);
set(gca,'FontSize',12,'Position',[0.06 0.6 0.43 0.4]);box on; grid on;
xlabel('(a) t = 0s','FontName','Times New Roman');
annotation('textbox',[.45 .61 .05 .05],'String','X','EdgeColor','none','FontSize',14);
annotation('textbox',[.06 .94 .05 .05],'String','Y','EdgeColor','none','FontSize',14)

subplot(2,2,2);
plot_agents(q_traj,qod,po,Nv,T/4,dT);
axis([-5.5   12.8   -7.5   8]);
set(gca,'FontSize',12,'Position',[0.56 0.6 0.43 0.4]);box on; grid on;
xlabel('(b) t = 20s','FontName','Times New Roman');
annotation('textbox',[.95 .61 .05 .05],'String','X','EdgeColor','none','FontSize',14);
annotation('textbox',[.56 .94 .05 .05],'String','Y','EdgeColor','none','FontSize',14)


subplot(2,2,3);
plot_agents(q_traj,qod,po,Nv,T/2,dT);box on; grid on;
axis([-5.5   12.8   -7.5   8]);
set(gca,'FontSize',12,'Position',[0.06 0.1 0.43 0.4]);
xlabel('(c) t = 40s','FontName','Times New Roman');
annotation('textbox',[.45 .11 .05 .05],'String','X','EdgeColor','none','FontSize',14);
annotation('textbox',[.06 .44 .05 .05],'String','Y','EdgeColor','none','FontSize',14)


subplot(2,2,4);
plot_agents(q_traj,qod,po,Nv,T,dT); box on; grid on;
axis([-5.5   12.8   -7.5   8]);
set(gca,'FontSize',12,'Position',[0.56 0.1 0.43 0.4]);
xlabel('(d) t = 80s','FontName','Times New Roman')
annotation('textbox',[.95 .11 .05 .05],'String','X','EdgeColor','none','FontSize',14);
annotation('textbox',[.56 .44 .05 .05],'String','Y','EdgeColor','none','FontSize',14);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf(figure(5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
    D = 1;
    for j = 1:N
        if i ~= j
            D = D.*dist(:,i,j);
        end
    end
    plot((1:Ns)*dT, D.^(1/(N-1)),'LineWidth',1); hold on; grid on;
end
set(gca,'FontSize',14,'Position',[0.09 0.2 0.88 0.77]);
xlabel('Time (s)','FontName','Times New Roman');
ylabel('$\Delta_i$','Interpreter','latex','Rotation',0);
axis([0 T 0 9]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf(figure(6)); i = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fd = axes('Position',[0.55 0.055 0.44 0.41]);
plot((1:Ns)*dT, up(1,:,i),'b','LineWidth',2); hold on; grid on;
plot((1:Ns)*dT, up(2,:,i),'r','LineWidth',2); hold on;
axis([0 T -0.6 0.8]); yticks([-0.5 0 0.5])
legend('$u_{1x}(t)$','$u_{1y}(t)$','Interpreter','latex','Location','best')
xlabel('(d)','FontName','Times New Roman');
set(fd,'FontSize',18);
% hold on;
% fd1 = axes('position',[0.38 0.17 0.4 0.07],'Box','on'); hold on
% plot((10/dT:30/dT)*dT,up(1,(10/dT:30/dT),i),'b','LineWidth',1); grid on;
% plot((10/dT:30/dT)*dT,up(2,(10/dT:30/dT),i),'r','LineWidth',1); grid on;
% axis([10 30 -1.1 1.5]); yticks([-1 -0.5 0 0.5 1.5]);
% set(fd1,'FontSize',10);


fc = axes('Position',[0.05 0.055 0.44 0.41]);
plot((1:Ns)*dT, uphi(:,i),'b','LineWidth',2); hold on; grid on;
axis([0 T -0.4 0.4]); yticks([-0.2 0 0.2])
legend('$u_{1\phi}(t)$','Interpreter','latex');
xlabel('(c)','FontName','Times New Roman');
set(fc,'FontSize',18);


% fc1 = axes('position',[0.45 0.345 0.4 0.07],'Box','on'); hold on
% plot((10/dT:30/dT)*dT,uphi((10/dT:30/dT),i),'b','LineWidth',1); grid on;
% axis([10 30 -0.4 0.41]);
% set(fc1,'FontSize',10);

fb = axes('Position',[0.55 0.55 0.44 0.41]);
plot((1:Ns)*dT, z_traj(1,1:Ns,i),'b','LineWidth',2); hold on; grid on;
plot((1:Ns)*dT, z_traj(2,1:Ns,i),'r','LineWidth',2);
plot((1:Ns)*dT, dqod(1,:) - (q_traj(2,1:Ns,i)- qod(2,:))*(-Ome),'m--','LineWidth',2);
plot((1:Ns)*dT, dqod(2,:) + (q_traj(1,1:Ns,i)- qod(1,:))*(-Ome),'c--','LineWidth',2);
axis([0 T -0.3 0.5])
legend('$v_{1x}(t)$','$v_{1y}(t)$',...
    '$\dot{x}_{r}(t)-(y_1(t)-y_{r}(t))\omega_1(t)$',...
    '$\dot{y}_{r}(t)+(x_1(t)-x_{r}(t))\omega_1(t)$',...
    'Interpreter','latex', 'NumColumns',2,...
    'Position',[0.360,0.670,0.574,0.062]);
xlabel('(b)','FontName','Times New Roman');
set(fb,'FontSize',18);


% subplot(4,1,1);
fa = axes('Position',[0.05 0.55 0.44 0.41]);
plot((1:Ns)*dT, z_traj(3,1:Ns,i),'b','LineWidth',2); hold on; grid on;
axis([0 T -0.12 0.01]);
yticks([-0.12 -0.08 -0.04 0]);
legend('$\omega_1(t)$','Interpreter','latex');
xlabel('(a)','FontName','Times New Roman');
set(fa,'FontSize',18);