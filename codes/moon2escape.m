% Copyright 2022 Nikita
% Created on: 15/10/22
% Outputs are all in the orbital frame of the moon, with x-axis pointing
% towards the ascending node of the moon

function [kepSC_flyby, nu_esc, flybyEt, kepM_flyby] = moon2escape(V_esc, kepMesc_F3, i_a, mu, r_esc, escapeEt, tol_OMEGA, tol_V)
%% MOON TO ESCAPE TRAJECTORY
v_esc = norm(V_esc);
if V_esc(3)>0
    i_z = 1;
else
    i_z = -1;
end
err = 1;
errV = 1;
convergence = 0;
iter = 1;
iterMax = 20;
[cartesianMesc_F3(1:3), cartesianMesc_F3(4:6)] = po2pv(kepMesc_F3,mu);
R_M = cartesianMesc_F3(1:3);

fprintf('\n The required escape velocity is: %fi + %fj + %fk', V_esc(1), V_esc(2), V_esc(3));
while convergence == 0 || iter <= iterMax
    if iter == iterMax
        cprintf('_Error','\nMaximum Iteration reached and no convergence reached.\n')
        fprintf('Convergence Report:\n')
        fprintf('Error in velocity: %f \n',errV);
        fprintf('Error in OMEGA_2: %f \n',err);
        break;
    end
    fprintf('\n======================================================================')
    fprintf("\nIteration: %d", iter);
    
    r_M = norm(R_M);
    U_n = i_a*R_M/r_M;
    OMEGA_2 = atan2(U_n(2),U_n(1));

    a_2 = -mu/(v_esc^2-2*mu/r_esc);
    alpha = acos(dot(U_n,V_esc)/v_esc);
    U_h = cross(U_n,V_esc)/norm(cross(U_n,V_esc));
    i_2 = acos(U_h(3));

    % get e2
    c_a = (a_2/r_M)^2;
    c_b = 2*(a_2/r_M)*(1-i_a*cos(alpha))-(sin(alpha))^2;
    c_c = (1-i_a*cos(alpha));
    if i_z*i_a == 1
        e_2 = max(sqrt(roots([c_a, c_b, c_c])+1));
    else
        e_2 = min(sqrt(roots([c_a, c_b, c_c])+1));
    end
    omega_2 = acos((a_2*(1-e_2.^2)/r_M-1)./e_2);
    PHI = acos(1./e_2);

    if i_a == +1
        nu_2 = -omega_2;
    else
        nu_2 = omega_2;
        omega_2 = -nu_2+pi;
    end
        
    if i_a*i_z == -1 && nu_2<-PHI
        cprintf('_Error',"\nThis flyby is not feasible.\n");
        break;
    end

    kepSC_F3 = [a_2, e_2, i_2, OMEGA_2, omega_2, nu_2];
    if i_a == 1
        nu_esc = acos((a_2*(1-e_2.^2)/r_esc-1)./e_2);
    else
        nu_esc = 2*pi-acos((a_2*(1-e_2.^2)/r_esc-1)./e_2);
    end

    TOFesc = trueAnomaly2time(kepSC_F3, nu_esc, mu)-trueAnomaly2time(kepSC_F3,nu_2, mu);   %rm to rsoi should be in days
    flybyEt = escapeEt-TOFesc;

    [cartesianSCesc_F3(1:3), cartesianSCesc_F3(4:6)] = po2pv([a_2, e_2, i_2, OMEGA_2, omega_2, nu_esc],mu);
    [cartesianSC2_F3(1:3), cartesianSC2_F3(4:6)] = po2pv([a_2, e_2, i_2, OMEGA_2, omega_2, nu_2],mu);
    [cartesianM2_F3(1:3), cartesianM2_F3(4:6)] = po2pv([kepMesc_F3(1:5); kepMesc_F3(6)-time2trueAnomaly(kepMesc_F3,TOFesc,mu)],mu);
    
    fprintf('\n The escape velocity of spacecraft is: %fi + %fj + %fk', cartesianSCesc_F3(4), cartesianSCesc_F3(5), cartesianSCesc_F3(6));
    fprintf('\n The spacecrafts location at flyby is: %fi + %fj + %fk', cartesianSC2_F3(1), cartesianSC2_F3(2), cartesianSC2_F3(3));
    fprintf('\n The moons location at flyby is      : %fi + %fj + %fk', cartesianM2_F3(1), cartesianM2_F3(2), cartesianM2_F3(3));
    fprintf('\n The moons location at escape is     : %fi + %fj + %fk', cartesianMesc_F3(1), cartesianMesc_F3(2), cartesianMesc_F3(3));
    kepM_F3 = pv2po(cartesianM2_F3(1:3)', cartesianM2_F3(4:6)', mu);
   
    R_M = cartesianM2_F3(1:3);
    r_M = norm(R_M);
    U_nNew = i_a*R_M/r_M;
    OMEGA_2New = atan2(U_nNew(2),U_nNew(1));

    err = abs(OMEGA_2New-OMEGA_2);
    errV = abs(norm(V_esc'-cartesianSCesc_F3(4:6))); 
    if errV < tol_V && err < tol_OMEGA
        fprintf("\n Convergence reached! \n");
        convergence = 1;
        break
    end
    OMEGA_2 = OMEGA_2New;
    iter = iter+1;
    
end
kepSC_flyby = kepSC_F3;
kepM_flyby = kepM_F3;
end