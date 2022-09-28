% Created on: 14/09/22 20:31

function [PO_Tf, PO_Ti, delV1, delV2] = getTransferOrbit(PO_A, PO_D, TOF, mu, Nrev, Ncase)

[RR_A,VV_A] = po2pv(PO_A,mu);
[RR_D,VV_D] = po2pv(PO_D,mu);
[A_T,P_T,E_T,ERROR_T,VI_T,VF_T,TPAR_T,THETA_T] = lambertMR(RR_D(1:3), RR_A(1:3), TOF, mu, 0, Nrev, Ncase,2);
%[VI_T, VF_T, extremal_distances, exitflag] = lambert(RR_D', RR_A', TOF, Nrev, mu);

delV1 = VI_T'-VV_D;
delV2 = VV_A-VF_T';

PO_Tf = pv2po(RR_A,VF_T', mu);
PO_Ti = pv2po(RR_D,VI_T', mu);
