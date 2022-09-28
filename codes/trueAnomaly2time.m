function t = trueAnomaly2time(OE, theta, mu)

e = OE(2);
a = OE(1);
if a<0
    a = -a;
end

if e>1
    type = 'hyperbola';
    h = sqrt(a*mu*(e^2-1));
elseif e == 0
    type = 'parabola';
    h = sqrt(a*mu*(1-e^2));
elseif e<1
    type = 'ellipse';
    h = sqrt(a*mu*(1-e^2));
elseif e == 1
    type = 'circle';
    h = sqrt(a*mu*(1-e^2));
end

switch type
    case 'hyperbola'
        F = log((sin(theta)*sqrt(e^2-1)+cos(theta)+e)/(1+e*cos(theta)));
        Mh = e*sinh(F)-F;
        t = (e^2-1)^(-3/2)*(h^3/mu^2)*Mh;
    case 'ellipse'
        E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
        Me = E - e*sin(E);
        t = (1-e^2)^(-3/2)*(h^3/mu^2)*Me;
    case 'circle'
        r = h^2/mu;
        T = 2*pi*r^(3/2)/sqrt(mu);
        t = theta/(2*pi)*T;
    case 'parabola'
        Mp = 1/2*tan(theta/2)+1/6*(tan(theta/2))^3;
        t = h^3/mu^2*Mp;
end

end