% Calculation of position and velocity vector in heliocentric coordinate
% system
function [rr, vv] = po2pv(PO, mu) 

% input: PO vector of classical orbital parameters 
%        mu gravity parameter in units consistent with a
% output: rr position vecotor in units consistent with mu and a
%         vv velocity vector in units consistent with mu and a

a = PO(1);
e = PO(2);
i = PO(3);
Om = PO(4);
om = PO(5);
theta = PO(6);

A = [cos(om+theta) -sin(om+theta) 0;
     sin(om+theta) cos(om+theta)  0;
                0             0   1];

if i<0
    i = pi+i;
end
    
B = [1      0       0;
     0 cos(i)  -sin(i);
     0 sin(i)   cos(i)];
 
C =  [cos(Om) -sin(Om) 0;
      sin(Om) cos(Om)  0;
          0       0   1];
      
p = a*(1-e^2);

r = [1/(1+e*cos(theta))*p 0 0]';
v = sqrt(mu/p)*[e*sin(theta) 1+e*cos(theta) 0]';

rr = C*B*A*r;
vv = C*B*A*v;