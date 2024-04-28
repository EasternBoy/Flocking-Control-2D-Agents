qi = [1;1;0];
qj = [2;2;0];
qod = [0;0;0];
poi = [0.4 0; 0 -0.2; -0.4 0; 0 0.2]';
poj = [0.4 0; 0 -0.2; -0.4 0; 0 0.2]';
Mi = 4;
Mj = 4;
dD = 1; 
dM = 5;
lam = 0.1;
mu = 100;

codegen -report GH_ij.m -args {qi,qj,qod,poi,poj,Mi,Mj,dD,dM,lam,mu}