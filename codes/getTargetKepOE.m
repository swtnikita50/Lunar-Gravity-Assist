
function [kep, cartesianSV] = getTargetKepOE(target, et, frame, observer, mu)


% Construct a meta kernel, "standard.tm”, which will be used to load the needed
% generic kernels: "naif0009.tls," "de432s.bsp," and "pck0009.tpc."
% Load the generic kernels using the meta kernel, and a Cassini spk.
%****************************%
%     Load SPICE kernels     %
%****************************%
%cspice_furnsh('./kernel.txt')
%****************************%

[cartesianSV,ltime] = cspice_spkezr(target, et, frame, 'NONE',observer);
kep = pv2po(cartesianSV(1:3), cartesianSV(4:6),mu);
end