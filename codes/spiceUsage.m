addpath('D:\NIKKY\Software\mice\lib')
addpath('D:\NIKKY\Software\mice\src\mice')

% Construct a meta kernel, "standard.tm”, which will be used to load the needed
% generic kernels: "naif0009.tls," "de421.bsp,” and "pck0009.tpc.”
% Load the generic kernels using the meta kernel, and a Cassini spk.
%****************************%
%     Load SPICE kernels     %
%****************************%
cspice_furnsh('./kernel.txt')
%****************************%

et = cspice_str2et('Mar 1, 2022');
[cartesianSV,ltime] = cspice_spkezr('MOON', times, 'J2000', 'NONE','EARTH');
%scatter3(pos(1,:), pos(2,:), pos(3,:));hold on;
%grid on
%scatter3(0,0,0,'p')
%muE = 3.986004418e5;    % km3/s2
% output1 = cspice_oscelt(cartesian,et,muE);
% cartesianSV = cartesianSV*1000; % convert to m from km
muE = 3.986004418e14;    % m3/s2
kep_J2000 = pv2po(cartesianSV(1:3), cartesianSV(4:6),muE);

% PQW = PQW2IJK(cartesianSV(1:3), cartesianSV(4:6), muE);
% cartesianSV_PQW(1:3) = PQW'*cartesianSV(1:3);
% cartesianSV_PQW(4:6) = PQW'*cartesianSV(4:6);
% cartesianSV_PQW = cartesianSV_PQW';