%% =========================================================================
%% 1. PARAMETER INITIALIZATION
%% =========================================================================
function QuadrotorPID()
% Physical parameters
m   = 1.121;         % Quadrotor mass (kg)
g   = 9.81;          % Gravity acceleration (m/s^2)

Ixx = 0.01;          % Moments of inertia (kg*m^2)
Iyy = 0.01;
Izz = 0.0148;
Ir  = 2.83e-5;       % Rotor moment of inertia

k   = 2.98e-5;       % Lift constant
b   = 3.23e-7;       % Drag constant
l   = 0.25;          % Arm length (m)

Ax  = 5.56e-4;       % Drag coefficients
Ay  = 5.56e-4;
Az  = 6.35e-4;

%% =========================================================================
%% 2. INITIAL CONDITIONS & SIMULATION TIMING
%% =========================================================================
% Initial state vector x0 (12 x 1)
% [x, y, z, vx, vy, vz, phi, theta, psi, phi_dot, theta_dot, psi_dot]

x0 = zeros(12, 1);  

t0 = 0;         
Tf = 20;        

%% =========================================================================
%% 3. DIFFERENTIAL EQUATION SOLVER (ODE45)
%% =========================================================================
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-5, 'MaxStep', 0.02);

% Function handle directly calls nested function 'f'
[t, X] = ode45(@f, [t0 Tf], x0, opts);

% Extract states
x = X(:,1);   y = X(:,2);   z = X(:,3);
phi = X(:,7); theta = X(:,8); psi = X(:,9);

%% =========================================================================
%% 4. REFERENCE TRAJECTORY GENERATION FOR PLOTTING
%% =========================================================================
xd = 0.5 * cos(pi * t / 5);
yd = 0.5 * sin(pi * t / 5);
zd = 1 - 0.5 * cos(pi * t / 5);

%% =========================================================================
%% 5. PLOTTING RESULTS
%% =========================================================================
figure(1);
plot3(x, y, z, 'b', 'LineWidth', 2); hold on;
plot3(xd, yd, zd, 'r--', 'LineWidth', 2);
grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
legend('Actual', 'Reference');
title('3D Trajectory Tracking');
axis equal; view(45, 30);

figure(2);
subplot(3,1,1); plot(t, x, 'b', t, xd, 'r--', 'LineWidth', 1.5); grid on;
ylabel('x (m)'); legend('actual','desired')
title('Position Tracking');
subplot(3,1,2); plot(t, y, 'b', t, yd, 'r--', 'LineWidth', 1.5); grid on;
ylabel('y (m)'); legend('actual','desired')
subplot(3,1,3); plot(t, z, 'b', t, zd, 'r--', 'LineWidth', 1.5); grid on;
ylabel('z (m)'); xlabel('Time (s)');  legend('actual','desired')

figure(3);
subplot(3,1,1); plot(t, rad2deg(phi), 'LineWidth', 1.5);
title('Roll (\phi)'); grid on; ylabel('deg');
subplot(3,1,2); plot(t, rad2deg(theta), 'LineWidth', 1.5);
title('Pitch (\theta)'); grid on; ylabel('deg');
subplot(3,1,3); plot(t, rad2deg(psi), 'LineWidth', 1.5);
title('Yaw (\psi)'); grid on; ylabel('deg'); xlabel('Time (s)');

%% =========================================================================
%% CONTROL FUNCTION
%% =========================================================================
function dX = f(t, X)
    persistent int_p int_a t_prev
    if isempty(int_p)
        int_p = zeros(3,1);
        int_a = zeros(3,1);
        t_prev = t;
    end
    
    dt = max(t - t_prev, 1e-6); % Dynamic time step
    t_prev = t;

    % Reference Trajectory
    xd_t = 0.5 * cos(pi * t / 5);
    yd_t = 0.5 * sin(pi * t / 5);
    zd_t = 1 - 0.5 * cos(pi * t / 5);
    ref = [xd_t; yd_t; zd_t];
    psi_d = 0;

    % State Unpacking
    pos   = X(1:3);
    vel   = X(4:6);
    ang   = X(7:9);
    rates = X(10:12);

    phi       = ang(1);
    theta     = ang(2);
    psi       = ang(3);
    phi_dot   = rates(1);
    theta_dot = rates(2);
    psi_dot   = rates(3);

    % External Disturbances
    dx     = 0.3 * sin(t);
    dy     = 0.3 * sin(t);
    dz     = 0.1 * cos(0.5 * t);
    dphi   = 0.2 * sin(2 * t);
    dtheta = 0.2 * sin(2 * t);
    dpsi   = 0.1 * sin(t);

    % PID Controller Parameters
    Kp_p = [6; 6; 7];
    Ki_p = [0.05; 0.05; 0.20];
    Kd_p = [4; 4; 7];

    Kp_a = [12; 12; 6];
    Ki_a = [0.02; 0.02; 0.01];
    Kd_a = [4; 4; 2];

%% =========================================================================
%% Outer Loop: Position Control
%% =========================================================================
    e_p = ref - pos;
    e_v = -vel; 
    int_p = int_p + e_p * dt;
    int_p = sat(int_p, 2);

    u = Kp_p .* e_p + Kd_p .* e_v + Ki_p .* int_p;
    ux = u(1); uy = u(2); uz = u(3);

    % Thrust and Target Roll/Pitch
    T = m * sqrt(ux^2 + uy^2 + (uz + g)^2); 
    
    phi_d   = asin((ux*sin(psi_d) - uy*cos(psi_d)) / sqrt(ux^2 + uy^2 + (uz + g)^2));

    theta_d = atan((ux*cos(psi_d) + uy*sin(psi_d)) / (uz + g));

    ang_ref = [phi_d; theta_d; psi_d];

%% =========================================================================
%% Inner Loop: Attitude Control
%% =========================================================================    
    e_a    = ang_ref - ang;
    e_rate = -rates;
    int_a  = int_a + e_a * dt;
    int_a  = sat(int_a, 1);

    tau       = Kp_a .* e_a + Kd_a .* e_rate + Ki_a .* int_a;
    tau_phi   = tau(1);
    tau_theta = tau(2);
    tau_psi   = tau(3);

    % Equations of Motion
    R = rotmat(phi, theta, psi);
    TB = [0; 0; T];
    G  = [0; 0; -m * g];
    FD = [Ax * vel(1); Ay * vel(2); Az * vel(3)];
    d_pos = [dx; dy; dz];

    acc = (R * TB + G - FD + d_pos) / m;

    % Motor speeds
    w1 = sqrt(max(0,T/(4*k)) - tau_theta/(2*k*l) - tau_psi/(4*b));
    w2 = sqrt(max(0,T/(4*k)) - tau_phi  /(2*k*l) + tau_psi/(4*b));
    w3 = sqrt(max(0,T/(4*k)) + tau_theta/(2*k*l) - tau_psi/(4*b));
    w4 = sqrt(max(0,T/(4*k)) + tau_phi  /(2*k*l) + tau_psi/(4*b));

    w_alpha = w1 - w2 + w3 - w4;

    % Angular Accelerations
    phi_dot2   = (tau_phi + (Iyy - Izz)*theta_dot*psi_dot + Ir*w_alpha*theta_dot + dphi)/Ixx;
    theta_dot2 = (tau_theta + (Izz - Ixx)*phi_dot*psi_dot - Ir*w_alpha*phi_dot + dtheta)/Iyy;
    psi_dot2   = (tau_psi + (Ixx - Iyy)*phi_dot*theta_dot + dpsi)/Izz;

    % Derivatives Vector Output
    dX = zeros(12,1);
    dX(1:3)  = vel;
    dX(4:6)  = acc;
    dX(7:9)  = rates;
    dX(10)   = phi_dot2;
    dX(11)   = theta_dot2;
    dX(12)   = psi_dot2;
end


%% =========================================================================
%% HELPER FUNCTION
%% =========================================================================
function R = rotmat(phi, theta, psi)
    cphi = cos(phi);     sphi = sin(phi);
    ctheta = cos(theta); stheta = sin(theta);
    cpsi = cos(psi);     spsi = sin(psi);

    R = [ctheta*cpsi,  sphi*stheta*cpsi - cphi*spsi,  cphi*stheta*cpsi + sphi*spsi;
         ctheta*spsi,  sphi*stheta*spsi + cphi*cpsi,  cphi*stheta*spsi - sphi*cpsi;
         -stheta,      sphi*ctheta,                   cphi*ctheta];
end

function y = sat(u, lim)
    y = min(max(u, -lim), lim);
end

end
