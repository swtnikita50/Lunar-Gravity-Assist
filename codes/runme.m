% Created on 20/09/2022

% tentative phi2, compute orbital parameters and TOF from moon to boundary
% of SOI

% Guess phi2
% Compute orbital parameters and TOF
% Actual position of moon at flyby
% get new phi2 from this position of the moon

%==========================================================================
%% TO BE CHECKED
% in MOON TO ESCAPE TRAJECTORY
% nu_2 <-PHI condition in the code, doesn't seem to work 
% complex OE parameters in descedning calculation for +ve Vescz
%==========================================================================

close all;
clear all;
clc;

addpath('D:\NIKKY\Software\mice\lib')
addpath('D:\NIKKY\Software\mice\src\mice')
cspice_furnsh('./kernel.txt')

muS = 1.32712440018e11;
muE = 3.986004418e5;    % km3/s2
muM = 7.34767309e22*6.6743e-20;
Nrev = 0; Ncase = 0;
r_SOI = 1e6;    %km


%% HELIOCENTRIC TRAJECTORY
escapeObject = 'EARTH';
arrivalObject = 'MARS BARYCENTER';
escapeEpoch = datestr(datetime('2001-01-04 10:05:24'));
escapeEt = cspice_str2et(escapeEpoch);
TOF = 210*60*60*24;
arrivalEt = escapeEt +TOF;
flybyEt = escapeEt;
arrivalEpoch = cspice_et2utc(arrivalEt, 'C', 1e-3);

fprintf("\n The Escape Time from %s is: %s", escapeObject ,escapeEpoch);
fprintf("\n The Arrival Time to %s is: %s", arrivalObject ,arrivalEpoch);
fprintf("\n The Time of Flight is %d day(s)", TOF/(60*60*24));
fprintf('\n======================================================================')
    
[kepDJ2000, cartesianDSV] = getTargetKepOE(escapeObject, escapeEt, 'J2000', 'SUN', muS);
[kepAJ2000, cartesianASV] = getTargetKepOE(arrivalObject, arrivalEt,'J2000', 'SUN', muS);
[PO_Tf, PO_Ti, delV1, delV2] = getTransferOrbit(kepAJ2000, kepDJ2000, TOF, muS, Nrev, Ncase);
plotOrbit(PO_Tf, PO_Ti(6), PO_Tf(6), muS, 100, "Transfer Trajectory", 'b');
plotOrbit(kepDJ2000, 0, 2*pi, muS, 100, "Initial Orbit", 'g');
plotOrbit(kepAJ2000, 0, 2*pi, muS, 100, "Target Orbit", 'r');
scatter3(0,0,0,696,'oy','filled','DisplayName','Center Body');
title('Earth Mars Transfer Trajectory in SunJ2000');
xlabel('x_{J2000} (km)'); ylabel('y_{J2000} (km)'); zlabel('z_{J2000} (km)');
legend;


%% MOON TO ESCAPE TRAJECTORY

tol_OMEGA = 1e-6;
tol_V = 1e-1;

[kepMescJ2000, cartesianMSVesc_F2] = getTargetKepOE('MOON', escapeEt, 'J2000', 'EARTH', muE);
ROT = [cos(-kepMescJ2000(5)), sin(-kepMescJ2000(5)), 0;...
    -sin(-kepMescJ2000(5)), cos(-kepMescJ2000(5)), 0;...
    0, 0, 1];
PQW = PQW2IJK(cartesianMSVesc_F2(1:3), cartesianMSVesc_F2(4:6), muS);
V_esc = ROT*PQW'*delV1; % convert delV1 to the PQW frame with omega_m rotation and you get V_esc;
pvTf = po2pv(PO_Tf,muS);
r_esc = r_SOI;

cartesianMSVesc_F3(1:3,1) = ROT*PQW'*cartesianMSVesc_F2(1:3);
cartesianMSVesc_F3(4:6,1) = ROT*PQW'*cartesianMSVesc_F2(4:6);
kepMesc_F3 = pv2po(cartesianMSVesc_F3(1:3), cartesianMSVesc_F3(4:6), muE);
[cartesianMesc_F3(1:3), cartesianMesc_F3(4:6)] = po2pv(kepMesc_F3,muE);
R_M = cartesianMesc_F3(1:3);

fprintf('\n ===================================================================== \n');
fprintf("\n Starting Analysis for Ascending node...\n");
fprintf('\n ===================================================================== \n');
[kepSC_F3, nu_esc, flybyEt, kepM_F3] = moon2escape(V_esc, kepMesc_F3, 1, muE, r_esc, escapeEt, tol_OMEGA, tol_V);
TOFesc = escapeEt-flybyEt;
figure(2)
title('Moon to Escape Trajectory: Ascending node flyby')
scatter3(0,0,0,20, 'ob',"filled", 'DisplayName', 'Earth');hold on;
grid on;
plotOrbit(kepMesc_F3, kepMesc_F3(6), kepMesc_F3(6)-time2trueAnomaly(kepMesc_F3,TOFesc,muE), muE, 100, "Moon's Orbit", 'r');
plotOrbit(kepSC_F3, kepSC_F3(6), nu_esc, muE, 100, "M2E Trajectory 1 Ascending", 'b'); hold on;
legend;

%% add perigee to moon

fprintf('\n ===================================================================== \n');
fprintf("\n Starting Analysis for Descending node...\n");
fprintf('\n ===================================================================== \n');
[kepSC_F3, nu_esc, flybyEt, kepM_F3] = moon2escape(V_esc, kepMesc_F3, -1, muE, r_esc, escapeEt, tol_OMEGA, tol_V);

