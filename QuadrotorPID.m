clear; clc; close all;

P.m   = 1.121;
P.g   = 9.81;

P.Ixx = 0.01;
P.Iyy = 0.01;
P.Izz = 0.0148;
P.Ir  = 2.83e-5;

P.k   = 2.98e-5;
P.b   = 3.23e-7;
P.l   = 0.25;

P.Ax  = 5.56e-4;
P.Ay  = 5.56e-4;
P.Az  = 6.35e-4;

P.Kp_pos = [6;6;7];
P.Ki_pos = [0.05;0.05;0.20];
P.Kd_pos = [4;4;7];

P.Kp_att = [12;12;6];
P.Ki_att = [0.02;0.02;0.01];
P.Kd_att = [4;4;2];
%%
X0 = zeros(12,1);   % [x y z vx vy vz phi theta psi phi_dot theta_dot psi_dot]
tspan = [0 40];

opts = odeset('RelTol',1e-3,'AbsTol',1e-5,'MaxStep',0.02);

[t,X] = ode45(@(t,X) quad_model(t,X,P), tspan, X0, opts);

x = X(:,1); y = X(:,2); z = X(:,3);
phi = X(:,7); theta = X(:,8); psi = X(:,9);

xd = 0.5*cos(pi*t/5);
yd = 0.5*sin(pi*t/5);
zd = 1-0.5*cos(pi*t/5);

figure;
subplot(3,1,1)
plot(t,x,'b',t,xd,'r--','LineWidth',1.5)
grid on; ylabel('x (m)')

subplot(3,1,2)
plot(t,y,'b',t,yd,'r--','LineWidth',1.5)
grid on; ylabel('y (m)')

subplot(3,1,3)
plot(t,z,'b',t,zd,'r--','LineWidth',1.5)
grid on; ylabel('z (m)')
xlabel('Time (s)')

figure;
subplot(3,1,1)
plot(t,rad2deg(phi),'LineWidth',1.5);
title('phi')
grid on
subplot(3,1,2)
plot(t,rad2deg(theta),'LineWidth',1.5)
title('theta')
grid on
subplot(3,1,3)
plot(t,rad2deg(psi),'LineWidth',1.5)
title('psi')
grid on
xlabel('Time (s)')
ylabel('deg')

figure;
plot3(x,y,z,'b','LineWidth',2); hold on
plot3(xd,yd,zd,'r--','LineWidth',2)
grid on
xlabel('x'); ylabel('y'); zlabel('z')
legend('Actual','Reference')
title('3D Trajectory')
axis equal
view(45,30)

function dX = quad_model(t,X,P)

persistent int_pos int_att

if isempty(int_pos)
    int_pos = zeros(3,1);
    int_att = zeros(3,1);
end

persistent t_prev
if isempty(t_prev)
    t_prev = t;
end

dt = 0.01;

x = X(1); y = X(2); z = X(3);
vx = X(4); vy = X(5); vz = X(6);

phi   = X(7);
theta = X(8);
psi   = X(9);

phi_dot = X(10);
theta_dot = X(11);
psi_dot = X(12);

pos   = [x; y; z];
vel   = [vx; vy; vz];
ang   = [phi; theta; psi];
rates = [phi_dot; theta_dot; psi_dot];
%%
xd = 0.5*cos(pi*t/5);
yd = 0.5*sin(pi*t/5);
zd = 1-0.5*cos(pi*t/5);

ref = [xd; yd; zd];
psi_d = 0;


dx = 0.3*sin(t);
dy = 0.3*sin(t);
dz = 0.1*cos(0.5*t);

dphi   = 0.2*sin(2*t);
dtheta = 0.2*sin(2*t);
dpsi   = 0.1*sin(t);

e_pos = ref - pos;
e_vel = -vel;

int_pos = int_pos + e_pos*dt;
int_pos = sat(int_pos,2);
%%
u = P.Kp_pos.*e_pos + P.Kd_pos.*e_vel + P.Ki_pos.*int_pos;

ux = u(1);
uy = u(2);
uz = u(3);
%%
T = P.m * sqrt(ux^2 + uy^2 + (uz + P.g)^2);

phi_d = asin((ux*sin(psi_d) - uy*cos(psi_d)) / sqrt(ux^2 + uy^2 + (uz + P.g)^2));

theta_d = atan((ux*cos(psi_d) + uy*sin(psi_d)) / (uz + P.g));

ang_ref = [phi_d; theta_d; psi_d];
%%
e_ang  = ang_ref - ang;
e_rate = - rates;

int_att = int_att + e_ang*dt;
int_att = sat(int_att,1);

tau = P.Kp_att.*e_ang + P.Kd_att.*e_rate + P.Ki_att.*int_att;

tau_phi   = tau(1);
tau_theta = tau(2);
tau_psi   = tau(3);
%%
R = rotmat(phi,theta,psi);

TB = [0;0;T];
G  = [0;0;- P.m*P.g];

FD = [P.Ax*vx;
      P.Ay*vy;
      P.Az*vz];

d_pos = [dx; dy; dz];

acc = (R*TB + G - FD + d_pos)/P.m;


w1 = sqrt(max(0,(T/(4*P.k)) - tau_theta/(2*P.k*P.l) - tau_psi/(4*P.b)));
w2 = sqrt(max(0,(T/(4*P.k)) - tau_phi /(2*P.k*P.l) + tau_psi/(4*P.b)));
w3 = sqrt(max(0,(T/(4*P.k)) + tau_theta/(2*P.k*P.l) - tau_psi/(4*P.b)));
w4 = sqrt(max(0,(T/(4*P.k)) + tau_phi /(2*P.k*P.l) + tau_psi/(4*P.b)));

w_alpha = w1 - w2 + w3 - w4;

phi_dot2 = (tau_phi + (P.Iyy-P.Izz)*theta_dot*psi_dot + P.Ir*theta_dot*w_alpha + dphi)/P.Ixx;

theta_dot2 = (tau_theta + (P.Izz-P.Ixx)*phi_dot*psi_dot - P.Ir*phi_dot*w_alpha + dtheta)/P.Iyy;

psi_dot2 = (tau_psi + (P.Ixx-P.Iyy)*phi_dot*theta_dot+ dpsi)/P.Izz;

%%
dX = zeros(12,1);

dX(1)=vx;
dX(2)=vy;
dX(3)=vz;

dX(4)=acc(1);
dX(5)=acc(2);
dX(6)=acc(3);

dX(7)=phi_dot;
dX(8)=theta_dot;
dX(9)=psi_dot;

dX(10)=phi_dot2;
dX(11)=theta_dot2;
dX(12)=psi_dot2;

end
