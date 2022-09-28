% Created on 20/09/2022

% tentative phi2, compute orbital parameters and TOF from moon to boundary
% of SOI

% Guess phi2
% Compute orbital parameters and TOF
% Actual position of moon at flyby
% get new phi2 from this position of the moon
close all;
clear all;
clc;

addpath('D:\NIKKY\Software\mice\lib')
addpath('D:\NIKKY\Software\mice\src\mice')
cspice_furnsh('./kernel.txt')

muS = 1.32712440018e11;
muE = 3.986004418e5;    % km3/s2
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
[kepMJ2000, cartesianMSV_F2] = getTargetKepOE('MOON', flybyEt, 'J2000', 'EARTH', muE);
OMEGA_2 = 2*807*pi/1000;

convergence = 0;
tol = 1e-6;
err = 1;
iter = 1;
while convergence == 0 && iter < 1000
    fprintf('\n======================================================================')
 
    fprintf("\n Iteration: %d", iter);
    ROT = [cos(-kepMJ2000(5)), sin(-kepMJ2000(5)), 0;...
        -sin(-kepMJ2000(5)), cos(-kepMJ2000(5)), 0;...
        0, 0, 1];
    PQW = PQW2IJK(cartesianMSV_F2(1:3), cartesianMSV_F2(4:6), muS);
    V_esc = ROT*PQW'*delV1; % convert delV1 to the PQW frame with omega_m rotation and you get V_esc;
    v_esc = norm(V_esc);
    pvTf = po2pv(PO_Tf,muS);
    r_esc = r_SOI;
    
    fprintf('\n The required escape velocity is: %fi + %fj + %fk', V_esc(1), V_esc(2), V_esc(3));

    cartesianMSV_F3(1:3,1) = ROT*PQW'*cartesianMSV_F2(1:3);
    cartesianMSV_F3(4:6,1) = ROT*PQW'*cartesianMSV_F2(4:6);
    kepM_F3 = pv2po(cartesianMSV_F3(1:3), cartesianMSV_F3(4:6), muE);

    r_M = norm(cartesianMSV_F2(1:3));
    
    if V_esc(3)>0
        isVesczPositive = 1;
        i_z = 1;
    else
        isVesczPositive = 0;
        i_z = -1;
    end
    
    a_2 = -muE/(v_esc^2-2*muE/r_esc);
    U_n = [cos(OMEGA_2); sin(OMEGA_2); 0];
    alpha = acos(dot(U_n,V_esc)/v_esc);
    U_h = cross(U_n,V_esc)/norm(cross(U_n,V_esc));
    i_2 = acos(U_h(3));

    % get e2
    i_a = [1,-1];
    for i=1:2
        c_a(i) = (a_2/r_M)^2;
        c_b(i) = 2*(a_2/r_M)*(1-i_a(i)*cos(alpha))-(sin(alpha))^2;
        c_c(i) = (1-i_a(i)*cos(alpha));

        e_2h(i) = sqrt(1+(-c_b(i)+sqrt(c_b(i)^2-4*c_a(i)*c_c(i)))/(2*c_a(i))); % higher solution for asymptote and escape on the same side
        nu_2h(i) = acos((a_2*(1-e_2h(i)^2)/r_M-1)/e_2h(i));
        PHIh(i) = acos(1/e_2h(i));

        e_2l(i) = sqrt(1+(-c_b(i)-sqrt(c_b(i)^2-4*c_a(i)*c_c(i)))/(2*c_a(i)));
        nu_2l(i) = acos((a_2*(1-e_2l(i)^2)/r_M-1)/e_2l(i));
        PHIl(i) = acos(1/e_2l(i));
    end

% 
%     if nu_2l>-PHIl && isVesczPositive
%         idx = 2;
%     elseif nu_2l>-PHIl && ~isVesczPositive
%         idx = 1;
%     end
    idx = 2;
    e_2 = e_2l(idx);
    nu_2 = nu_2l(idx);
    PHI = PHIl(idx);
    
    if isVesczPositive
        %omega_2 = alpha+PHI-pi;
        omega_2 = -nu_2;
    else
        %omega_2 = pi-alpha+PHI;
        omega_2 = -nu_2+pi;
    end
    kepSC_F3 = [a_2, e_2, i_2, OMEGA_2, omega_2, nu_2];

    nu_esc = acos((a_2*(1-e_2^2)/r_esc-1)/e_2);
    TOFesc = trueAnomaly2time(kepSC_F3, nu_esc, muE)-trueAnomaly2time(kepSC_F3,nu_2, muE);%rm to rsoi should be in days
    flybyEt = escapeEt-TOFesc;
    
    [cartesianSCesc_F3(1:3), cartesianSCesc_F3(4:6)] = po2pv([a_2, e_2, i_2, OMEGA_2, omega_2, nu_esc],muE);
    [cartesianSC2_F3(1:3), cartesianSC2_F3(4:6)] = po2pv([a_2, e_2, i_2, OMEGA_2, omega_2, nu_2],muE);
    %[cartesianMesc_F3(1:3), cartesianMesc_F3(4:6)] = po2pv([kepM_F3(1:5); kepM_F3(6)+time2trueAnomaly(kepM_F3,TOFesc,muE)],muE);
    [cartesianMesc_F3(1:3), cartesianMesc_F3(4:6)] = po2pv([kepM_F3(1:5); kepM_F3(6)],muE);
    [cartesianM2_F3(1:3), cartesianM2_F3(4:6)] = po2pv(kepM_F3,muE);
    
    fprintf('\n The escape velocity of spacecraft at SOI is: %fi + %fj + %fk', cartesianSCesc_F3(4), cartesianSCesc_F3(5), cartesianSCesc_F3(6));
    fprintf('\n The flyby position of spacecraft is: %fi + %fj + %fk', cartesianSC2_F3(1), cartesianSC2_F3(2), cartesianSC2_F3(3));
    fprintf('\n The flyby position of moon is: %fi + %fj + %fk', cartesianM2_F3(1), cartesianM2_F3(2), cartesianM2_F3(3));

    fprintf('\n The time of fight from MOON to ESCAPE is: %f sec(s)', TOFesc);
    kepMJ2000 = getTargetKepOE('MOON', flybyEt, 'J2000', 'EARTH', muE);
    %OMEGA_2N = kepMJ2000(5)+kepMJ2000(6);%+kepMJ2000(3);
    %OMEGA_2N = -OMEGA_2N;
    err(iter) = abs(norm(V_esc'-cartesianSCesc_F3(4:6))) + abs(norm(cartesianSC2_F3(1:3)-cartesianM2_F3(1:3)));
    %err = abs(OMEGA_2N-OMEGA_2);
    if err < tol
        convergence = 1;
        break
    end
    %OMEGA_2 = OMEGA_2N;
    OMEGA_2 = OMEGA_2+0.0001;
    iter = iter+1;
end


gamma_esc = atan(e_2*sin(nu_esc/(1+e_2*cos(nu_esc))));
beta = pi/2-nu_esc+gamma_esc-PHI;

figure(2)
plotOrbit(kepSC_F3, nu_2, nu_esc, muE, 100, "Moon to Escape Trajectory", 'k'); hold on;
scatter3(0,0,0,20, 'ob',"filled", 'DisplayName', 'Earth');
grid on;
plotOrbit(kepM_F3, 0, 2*pi, muE, 100, "Moon's Orbit", 'r')
plotOrbit(kepM_F3, kepM_F3(6), kepM_F3(6)+ time2trueAnomaly(kepM_F3,TOFesc,muE), muE, 100, "Moon's Orbit", 'b')
[x,y,z]=sphere; x = r_SOI*x;y = r_SOI*y;z = r_SOI*z;
surf(x,y,z,'FaceAlpha',.2,'EdgeColor','none', "faceColor", 'g', 'displayName', "Earth's Sphere of Influence");
legend;
fprintf('\n ===================================================================== \n');

figure(3)
plotOrbit(kepSC_F3, 0, 2*pi, muE, 100, "Moon to Escape Trajectory", 'k'); hold on;
plotOrbit(kepM_F3, 0, 2*pi, muE, 100, "Moon's Orbit", 'r')
plotOrbit(kepM_F3, kepM_F3(6), kepM_F3(6)+ time2trueAnomaly(kepM_F3,TOFesc,muE), muE, 100, "Moon's Orbit", 'b')

figure(4)
plot(err);


%%
    i_a = [1,-1];
    for i=1:2
        c_a(i) = (a_2/r_M)^2;
        c_b(i) = 2*(a_2/r_M)*(1-i_a(i)*cos(alpha))-(sin(alpha))^2;
        c_c(i) = (1-i_a(i)*cos(alpha));

        e_2h(i) = sqrt(1+(-c_b(i)+sqrt(c_b(i)^2-4*c_a(i)*c_c(i)))/(2*c_a(i))); % higher solution for asymptote and escape on the same side
        nu_2h(i) = acos((a_2*(1-e_2h(i)^2)/r_M-1)/e_2h(i));
        PHIh(i) = acos(1/e_2h(i));

        e_2l(i) = sqrt(1+(-c_b(i)-sqrt(c_b(i)^2-4*c_a(i)*c_c(i)))/(2*c_a(i)));
        nu_2l(i) = acos((a_2*(1-e_2l(i)^2)/r_M-1)/e_2l(i));
        PHIl(i) = acos(1/e_2l(i));
    end