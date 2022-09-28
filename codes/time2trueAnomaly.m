function theta = time2trueAnomaly(OE, t, mu)

e = OE(2);
a = OE(1);

if e > 1
    type = 'hyperbola';
    h = sqrt(a*mu*(e^2-1));
elseif e == 0
    type = 'circle';
    h = sqrt(a*mu*(1-e^2));
elseif e < 1
    type = 'ellipse';
    h = sqrt(a*mu*(1-e^2));
elseif e == 1
    type = 'parabola';
    h = sqrt(a*mu*(1-e^2));
end

switch type
    case 'hyperbola'
        T = 2*pi*h^3/(mu^2*(e^2-1));
        Mh = 2*pi*t/T;
        theta = meanAnomaly2trueAnomaly(Mh, e);
    case 'ellipse'
        T = 2*pi*h^3/(mu^2*(1-e^2));
        Me = 2*pi*t/T;
        theta = meanAnomaly2trueAnomaly(Me, e);
    case 'circle'
        r = h^2/mu;
        T = 2*pi*r^(3/2)/sqrt(mu);
        theta = 2*pi*t/T;
    case 'parabola'
        T = mu^2/h^3;
        Mp = 2*pi*t/T;
        theta = meanAnomaly2trueAnomaly(Mp, e);
end

end
