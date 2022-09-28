function theta = meanAnomaly2trueAnomaly(M, e)

if e > 1
    type = 'hyperbola';
elseif e == 1
    type = 'parabola';
elseif e < 1
    type = 'ellipse';
end

tol = 1e-3;

switch type
    case 'hyperbola'
        F = 0;
        err = 1;
        while err> tol
            E = e*sin(F)-F-M;
            dEdF = e*cos(F)-F;
            Fnew = F-E/dEdF;
            err = abs(Fnew-F);
            F = Fnew;
        end
        theta = acos((cos(F)-e)/(1-e*cos(F)));
    case 'ellipse'
        E = 0;
        err = 1;
        while err> tol
            F = E-e*sin(E)-M;
            dFdE = 1-e*cos(E);
            Enew = E-F/dFdE;
            err = abs(Enew-E);
            E = Enew;
        end
        theta = acos((e-cos(E))/(e*cos(E)-1));
    case 'parabola'
        theta = 2*atan((3*M+sqrt(9*M^2+1))^(1/3) - (3*M+sqrt(9*M^2+1))^(-1/3));

end
