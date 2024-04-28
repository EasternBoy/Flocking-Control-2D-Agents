function plot_agents(q_traj,qod,po,Nv,T,dT)

hold on;
axis equal;

N = size(q_traj)*[0;0;1];


% Circle plot
r = 0.1;
xc1 = 0.3; yc1 = 0;

m = round(T/dT);

qT = zeros(3,N);

for k = 1:N
    qT(:,k) = q_traj(:,m,k);     % Position at T
end
    
for j = 1:12
    plot(q_traj(1,1:m,j),q_traj(2,1:m,j),'b.','LineWidth',1,'MarkerSize',3);
end
    
for j = 13:20
    plot(q_traj(1,1:m,j),q_traj(2,1:m,j),'m.','LineWidth',1,'MarkerSize',3);
end
    
for j = 21:N
    plot(q_traj(1,1:m,j),q_traj(2,1:m,j),'g.','LineWidth',1,'MarkerSize',3);
end

plot(qod(1,1:m), qod(2,1:m),'k-.','LineWidth',2); 

for j = 1:12
    phi = qT(3,j);
    poi = po(:,:,j)*[eye(Nv(j));zeros(4-Nv(j),Nv(j))];
    x = [cos(phi) -sin(phi)]*[poi poi] + qT(1,j)*ones(1,2*Nv(j));
    y = [sin(phi)  cos(phi)]*[poi poi] + qT(2,j)*ones(1,2*Nv(j));
    fill(x,y,'b','LineWidth',0.5);  
    xc = [cos(phi) -sin(phi)]*[xc1;yc1] + qT(1,j);
    yc = [sin(phi)  cos(phi)]*[xc1;yc1] + qT(2,j);
    circle(xc,yc,r);
end

for j = 13:20
    phi = qT(3,j);
    poi = po(:,:,j)*[eye(Nv(j));zeros(4-Nv(j),Nv(j))];
    x = [cos(phi) -sin(phi)]*[poi poi] + qT(1,j)*ones(1,2*Nv(j));
    y = [sin(phi)  cos(phi)]*[poi poi] + qT(2,j)*ones(1,2*Nv(j));
    fill(x,y,'m','LineWidth',0.5);  
    xc = [cos(phi) -sin(phi)]*[xc1;yc1] + qT(1,j);
    yc = [sin(phi)  cos(phi)]*[xc1;yc1] + qT(2,j);
    circle(xc,yc,r);
end

for j = 21:N
    phi = qT(3,j);
    poi = po(:,:,j)*[eye(Nv(j));zeros(4-Nv(j),Nv(j))];
    x = [cos(phi) -sin(phi)]*[poi poi] + qT(1,j)*ones(1,2*Nv(j));
    y = [sin(phi)  cos(phi)]*[poi poi] + qT(2,j)*ones(1,2*Nv(j));
    fill(x,y,'g','LineWidth',0.5);  
    xc = [cos(phi) -sin(phi)]*[xc1;yc1] + qT(1,j);
    yc = [sin(phi)  cos(phi)]*[xc1;yc1] + qT(2,j);
    circle(xc,yc,r);
end

hold off;
    