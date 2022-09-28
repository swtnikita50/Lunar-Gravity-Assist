% PQW to IJK frame
% Created on 23/09/22 01:09
% output is a rotation vector which rotates from PQW (perifocal) to IJK (inertial) frame
function PQW = PQW2IJK(pos_IJK, vel_IJK, mu)

e = (((norm(vel_IJK))^2-mu/norm(pos_IJK))*pos_IJK...
    -dot(pos_IJK,vel_IJK)*vel_IJK)/mu;

P = e/norm(e);
W = cross(pos_IJK,vel_IJK)/norm(cross(pos_IJK,vel_IJK));
Q = cross(W,P);

PQW = [P, Q, W];
end
