clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T  = 80;
dT = 0.4e-2;
Ns = T/dT;

N = 25;          % Number of agents
q_traj = zeros(3,Ns+1,N);
z_traj = zeros(3,Ns+1,N);
up     = zeros(2,Ns,N);
uphi   = zeros(Ns,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_ip   = 0.1;
k_iv   = 0.5;
k_iphi = 0.1;
k_iom  = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Potential functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dD  = 0.8;
dM  = 5;         % Communication range R = 3;
lam = 1e-3;
mu  = 100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0 = 4;
L1 = 12;
L2 = 8;
L3 = N - L1 - L2;
for i=1:1:L1
    q_traj(1,1,i) = R0*cos(2*pi*(i-1)/L1);
    q_traj(2,1,i) = R0*sin(2*pi*(i-1)/L1);
%     q_traj(3,1,i) = pi*(L1-i)/L1;
    q_traj(3,1,i)  = rand(1)*pi;
end
for i=1:1:L2
    q_traj(1,1,i+L1) = R0*2/3*cos(2*pi*(i-1)/L2);
    q_traj(2,1,i+L1) = R0*2/3*sin(2*pi*(i-1)/L2);
%     q_traj(3,1,i+L1) = pi*i/L2;
    q_traj(3,1,i+L1) = rand(1)*pi; 
end
for i=1:1:L3
    q_traj(1,1,i+L1+L2) = R0*1/3*cos(2*pi*(i-1)/L3);
    q_traj(2,1,i+L1+L2) = R0*1/3*sin(2*pi*(i-1)/L3);
%     q_traj(3,1,i+L1+L2) = pi*(L2-i)/L3;
    q_traj(3,1,i+L1+L2) = rand(1)*pi; 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Agent dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nv = 4*ones(1,N);  % 4-vertices agents
po = zeros(2,4,N);
for i = 1:N
    po(:,:,i)  = [0.4 0; 0 -0.2; -0.4 0; 0 0.2]';
end
S    = [0 -1;1 0];
dist = zeros(Ns,N,N);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:1:Ns
             
    q = zeros(3,N);
    z = zeros(3,N);
    for i = 1:N
        q(:,i) = q_traj(:,k,i);
        z(:,i) = z_traj(:,k,i);
    end

    M    = zeros(N,2);
    L    = zeros(N,1);

    for i = 1:N
        for j = i+1:1:N
            [d, G, H] = GH_ij_mex(q(:,i),q(:,j),qod(:,k),po(:,:,i),po(:,:,j),Nv(i),Nv(j),dD,dM,lam,mu);
            M(i,:) = M(i,:) + G;
            L(i)   = L(i)   + H;
            dist(k,i,j) = d;
        end

        for j = 1:1:i-1
            [d, G, H] = GH_ij_mex(q(:,j),q(:,i),qod(:,k),po(:,:,j),po(:,:,i),Nv(j),Nv(i),dD,dM,lam,mu);
            M(i,:) = M(i,:) - G;
            L(i)   = L(i)   - H;
            dist(k,i,j) = d;
        end   
        
        
        x  = q(1,i); y  = q(2,i); phi  = q(3,i);
        vx = z(1,i); vy = z(2,i); om   = z(3,i);
        
        xod  = qod(1,k);  yod  = qod(2,k); 
        dxod = dqod(1,k); dyod = dqod(2,k);
        ddxod = ddqod(1,k); ddyod = ddqod(2,k);
        
        alx = dxod - (y - yod)*om;
        aly = dyod + (x - xod)*om;
        
        ep1 = 0.2;
        ep2 = 1.2;

        uphi(k,i) = - L(i) + ddqod(3,k) - k_iphi*sat(phi-qod(3,k),ep1) - k_iom*sat(om - dqod(3,k),ep1);
        
        up(1,k,i) = - M(i,1) + ddxod - k_ip*(x - xod)/sqrt((x-xod)^2 + (y-yod)^2 + ep2^2)...
            - (vy-dyod)*om - (y-yod)*uphi(k,i) - k_iv*sat(vx-alx,1);
        
        up(2,k,i) = - M(i,2) + ddyod - k_ip*(y - yod)/sqrt((x-xod)^2 + (y-yod)^2 + ep2^2)...
            + (vx-dxod)*om + (x-xod)*uphi(k,i) - k_iv*sat(vy-aly,1);
    end

    for i = 1:N
        qz = [q(:,i); z(:,i)] + ([z(:,i); 0;0;0] + [0;0;0;up(:,k,i);uphi(k,i)])*dT;
        q_traj(:,k+1,i) = qz(1:3);
        z_traj(:,k+1,i) = qz(4:6);
    end

end

Figures


