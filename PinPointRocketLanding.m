%% Gagandeep Thapar 
% AERO 560: Final
% Convex Optimization Implementation:
        % Convex Optimization for Trajectory Generation: A Tutorial on Generating
        % Dynamically Feasible Trajectories Reliably and Efficiently

%% Housekeeping

clc;
close all;
clear all;

%% Setup

% setup world information (Earth)
world.x_hat = [1,0,0];
world.y_hat = [0,1,0];
world.z_hat = [0,0,1];

world.g = -9.807;
world.T = 23*3600 + 56*60 + 4.091;  % [s] sidereal time of earth
world.lat = 34.7420;
world.t_step = 1;

% create rocket struct (modeled off of Falcon 9/Merlin Engine for Earth landing)
rocket.g = -3.7114*world.z_hat;
rocket.omega = (2*pi/world.T)*(world.x_hat*cosd(world.lat) + world.y_hat*0 + world.z_hat*sind(world.lat));
rocket.m_dry = 1505;
rocket.m_wet = 1905;
rocket.isp = 225;
rocket.phi = 27; % [deg] max cant angle
rocket.alpha = 1/(rocket.isp*9.807*cosd(rocket.phi)); 
rocket.T_max = 6 * 3100;
rocket.thrust_min = 0.3 * rocket.T_max * cosd(rocket.phi);
rocket.thrust_max = 0.8 * rocket.T_max * cosd(rocket.phi);
rocket.approach_max = 86;  % [deg]
rocket.pointing_max = 40;  % [deg]
rocket.vel_max = 500 * 1000/3600;   % [m/s]
rocket.pos_initial = (2*world.x_hat' + 0*world.y_hat' + 1.5*world.z_hat')*1000;
rocket.vel_initial = 80*world.x_hat' + 30*world.y_hat' - 75*world.z_hat';

rocket.w_x = skew(rocket.omega);
rocket.A_c = [zeros(3,3), eye(3), zeros(3,1);
                -(rocket.w_x)^2, -2*rocket.w_x, zeros(3,1);
                zeros(1,7)];
rocket.B_c = [zeros(3,4);
                eye(3), zeros(3,1);
                zeros(1,3), -rocket.alpha];
rocket.p_c = [zeros(3,1);rocket.g';0];
[rocket.n, rocket.m] = size(rocket.B_c);
tFlight = (rocket.m_wet - rocket.m_dry)/(rocket.alpha*rocket.thrust_min);

%% PDG FFT

delt = 1;   % [sec] step time
tFlight = round(tFlight);
t = 0:delt:tFlight;

simSz = (tFlight)/delt + 1;

r = zeros(simSz,3);
v = zeros(simSz,3);
m = zeros(simSz,1);

r(1,:) = rocket.pos_initial;
v(1,:) = rocket.vel_initial;
m(1) = rocket.m_wet;

% scaling variables
s_r = zeros(1,3);
S_r = diag([max(1, abs(r(1,1))), max(1, abs(r(1,2))), max(1, abs(r(1,3)))]);
s_v = zeros(1,3);
S_v = diag([max(1, abs(v(1,1))), max(1, abs(v(1,2))), max(1, abs(v(1,3)))]);
s_z = (log(rocket.m_dry) + log(rocket.m_wet))/2;
S_z = log(rocket.m_wet)-s_z;
s_u = [0,0,0.5*(rocket.thrust_min/rocket.m_wet*cosd(rocket.pointing_max)+rocket.thrust_max/rocket.m_dry)];
S_u = diag([rocket.thrust_max/rocket.m_dry*sind(rocket.pointing_max), ...
            rocket.thrust_max/rocket.m_dry*sind(rocket.pointing_max), ...
            rocket.thrust_max/rocket.m_dry-s_u(3)]);
s_zeta = s_u(3);
S_zeta = S_u(3,3);

% glideslope information
H_gs = [cosd(rocket.approach_max), 0, -sind(rocket.approach_max);
        -cosd(rocket.approach_max), 0, -sind(rocket.approach_max);
        0, cosd(rocket.approach_max), -sind(rocket.approach_max);
        0, -cosd(rocket.approach_max), -sind(rocket.approach_max)];

h_gs = zeros(4,1);

% PDG setup
pdg = optimproblem;
r_s = optimvar('r_s', 3, simSz);
v_s = optimvar('v_s', 3, simSz);
z_s = optimvar('z_s', 1, simSz);
u_s = optimvar('u_s', 3, simSz-1);
zeta_s = optimvar('zeta_s', 1, simSz-1);

% Scaling
r = (S_r*r_s) + repmat(s_r, simSz,1)';
v = (S_v*v_s) + repmat(s_v, simSz,1)';
z = (S_z*z_s) + repmat(s_z, simSz,1)';
u = (S_u*u_s) + repmat(s_u, simSz-1,1)';
zeta = (S_zeta*zeta_s) + repmat(s_zeta, simSz-1,1)';

% initial and final constraints
r0 = optimconstr(1);
v0 = optimconstr(2);
z0 = optimconstr(3);

rF = optimconstr(4);
vF = optimconstr(5);
zF = optimconstr(6);

r0 = r(:,1) == rocket.pos_initial;
rF = r(:,end) == zeros(3,1);
v0 = v(:,1) == rocket.vel_initial;
vF = v(:,end) == zeros(3,1);
z0 = z(1) == log(rocket.m_wet);
zF = z(end) >= log(rocket.m_dry);

pdg.Constraints.r0 = r0;
pdg.Constraints.v0 = v0;
pdg.Constraints.z0 = z0;
pdg.Constraints.rF = rF;
pdg.Constraints.vF = vF;
pdg.Constraints.zF = zF;

dynamicConstraint = optimconstr(simSz-1,7);
thrustMinConstraint = optimconstr(simSz-1, 1);
thrustMaxConstraint = optimconstr(simSz-1, 1);
massMinConstraint = optimconstr(simSz, 1);
massMaxConstraint = optimconstr(simSz, 1);
glideslopeConstraint = optimconstr(simSz, 4);
attitudeConstraint = optimconstr(simSz-1,3);

[A,B,p] = discretize(rocket, delt);

% time-based constraints
for i = 1:simSz-1

    % dynamic constraints
    dynamicConstraint(i,1:3) = r(:,i+1)    == r(:,i) + (v(:,i))*delt;
    dynamicConstraint(i,4:6) = v(:,i+1)    == v(:,i) + ((zeta(i)*u(:,i)/z(i)) + rocket.g')*delt;
    dynamicConstraint(i,7)   = z(i+1)      == (z(i)) - log(rocket.alpha*zeta(i)*delt);

    % thrust
    zNOT = log(rocket.m_wet - rocket.alpha*rocket.thrust_max*t(i));
    uMin = rocket.thrust_min*exp(-zNOT);
    uMax = rocket.thrust_max*exp(-zNOT);
    sigZ = z(i) - zNOT;

    thrustMinConstraint(i) = zeta(i) >= uMin*(1-sigZ + 0.5*sigZ^2);
    thrustMaxConstraint(i) = zeta(i) <= uMax*(1-sigZ);

    % mass
    massMinConstraint(i) = zNOT <= z(i);
    massMaxConstraint(i) = z(i) <= log(rocket.m_wet - rocket.alpha*rocket.thrust_min*t(i));

    % glideslope
    glideslopeConstraint(i,:) = H_gs*(r(:,i)) <= h_gs;

    % attitude
    attitudeConstraint(i,:) = dot(u(:,i), [0;0;1]) >= zeta(i)*cosd(rocket.pointing_max);

end

% final state constraints
% mass
    zNOT = log(rocket.m_wet - rocket.alpha*rocket.thrust_max*t(simSz));
    massMinConstraint(end) = zNOT <= z(end);
    massMaxConstraint(end) = z(end) <= log(rocket.m_wet - rocket.alpha*rocket.thrust_min*t(end));

% glideslope
    glideslopeConstraint(end,:) = H_gs*(r(:,end)) <= h_gs;

disp('post')

pdg.Constraints.dynamicConstraint = dynamicConstraint;
pdg.Constraints.thrustMinConstraint = thrustMinConstraint;
pdg.Constraints.thrustMaxConstraint = thrustMaxConstraint;
pdg.Constraints.massMinConstraint = massMinConstraint;
pdg.Constraints.massMaxConstraint = massMaxConstraint;
pdg.Constraints.attitudeConstraint = attitudeConstraint;
pdg.Constraints.glideslopeConstraint = glideslopeConstraint;

pdg.Objective = sum(zeta);

%%

x0.r_s = zeros(3, simSz);
x0.u_s = zeros(3, simSz-1);
x0.v_s = zeros(3, simSz);
x0.z_s = zeros(1, simSz);
x0.zeta_s = zeros(1, simSz-1);

x0.r_s(:,1) = ones(3,1);
x0.v_s(:,1) = [1;1;-1];
x0.u_s(:,1) = zeros(3,1);
x0.z_s(1) = log(rocket.m_wet);
x0.zeta_s(1) = 0;

%%
solveFlag = 0;

tic
if solveFlag == 1

    disp('Solving problem with fmincon')
    
    options = optimoptions('fmincon', 'ConstraintTolerance',1e-3, 'OptimalityTolerance',1e-3, 'EnableFeasibilityMode',true, 'SubproblemAlgorithm','cg', 'Display', 'off');
    x = solve(pdg, x0, solver='fmincon');

else

    fprintf('Loading in data from mat file\n\n');
    load('subsurface_data.mat');

end
toc


%% Plots
close all

scaledr = S_r*x.r_s;

figure
subplot(3,1,1)
plot(t, scaledr(1,:))
ylabel('X [m]')

subplot(3,1,2)
plot(t, scaledr(2,:))
ylabel('Y [m]')

subplot(3,1,3)
plot(t, scaledr(3,:))
ylabel('Z [m]')
xlabel('Time [sec]')

sgtitle('Landing Trajectory')

figure
hold on
plot(scaledr(1,:),scaledr(3,:))
scatter(scaledr(1,1), scaledr(3,1), 'rx')
scatter(scaledr(1, end), scaledr(3,end), 'kx')
hold off
title('Downrange Trajectory')
legend('Trajectory', 'Starting Point', 'Landing Point')
xlabel('Distance from Landing Point [m]')
ylabel('Altitude [m]')

figure
hold on
plot(scaledr(2,:), scaledr(3,:))
scatter(scaledr(2,1), scaledr(3,1), 'rx')
scatter(scaledr(2, end), scaledr(3,end), 'kx')
hold off
title('Crossrange Trajectory')
xlabel('Crossrange [m]')
ylabel('Altitude [m]')
legend('Trajectory', 'Starting Point', 'Landing Point')

figure
plot3(scaledr(1,:), scaledr(2,:), scaledr(3,:))
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('3D Landing Trajectory')

vnorm = vecnorm(S_v*x.v_s,2,1);

figure
plot(t, vnorm)
title('Rocket Norm of Velocity')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')


mass =  rocket.m_dry*exp(S_z*x.z_s);
figure
hold on
plot([t(1),t(end)], [rocket.m_dry, rocket.m_dry], 'r--')
plot(t, mass, 'b')
hold off
title('Mass of Rocket')
legend('Rocket Dry Mass')
xlabel('Time [sec]')
ylabel('Mass [kg]')

%% Functions

function wx = skew(w)

    wx = [0, -w(3) w(2);
            w(3) 0 -w(1);
            -w(2) w(1) 0];

end

function [A,B,p] = discretize(rocket, delT)

    Ar = rocket.A_c;
    Br = rocket.B_c;
    pr = rocket.p_c;
    
    matr = [Ar, Br, pr;
            zeros(rocket.m+1, rocket.n+rocket.m+1)] * delT;
    
    M = exp(matr);

    A = M(1:rocket.n, 1:rocket.n);
    B = M(1:rocket.n, rocket.n+1:rocket.m+rocket.n);
    p = M(1:rocket.n, rocket.n+rocket.m+1);

end

