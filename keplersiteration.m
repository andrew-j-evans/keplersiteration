clear;
clc;
format longG;

dtv  = [0.1;1;10;60;300]; %Different step sizes used throughout the iteration
T = 1; %(days)

%Initial position and velocity vectors. Can be set to anything.  
R0 = [1131.34; -2282.343; 6672.423];
V0 = [-5.64305; 4.30333; 2.42879];

[a,e,i,RA,w,nu,mu] = RV2COE(R0,V0);

n = sqrt(mu / a^3);
%E = acos(deg2rad((e+cosd(nu))/(1+e*cosd(nu))));
%M = E - e*sin(E);

M = n * T*60*60*24;

if -pi<M<0
    E0=M-e;
elseif M>pi
    E0=M-e;
elseif 0<M<pi
    E0=M+e;
else
    E0=M(j);
end

error = 10;

while error > 10e-6
    E=E0+(M-E0+e.*sin(E0))./(1-e.*cos(E0));
    error = max(max(abs(E-E0)));
    E0=E;
end

nuF = acos((cos(E) - e)/(1-e*cos(E)));
nuFDeg = rad2deg(nuF);

[KRF, KVF] = COE2RV(a,e,i,RA,w,nuFDeg,mu);

KSME = norm(KVF)^2/2-mu/norm(KRF);
APK = cross(KRF, KVF);

%% Tables and Outputs
format shortG;
dt = table(reshape(dtv',1,[]),'VariableNames',"dt");
disp(dt)



finalRKep = table(KRF,'VariableNames',"Kepler's Final R Vector (km)");
finalVKep = table(KVF,'VariableNames',"Kepler's Final V Vector (km/s)");
disp(finalRKep);
disp(finalVKep);

orbitalElements = table(a, e, i, RA, w, nu, 'VariableNames',["a (km)","e","i (deg)","RA (deg)","w (deg)","nu (deg)"]);
disp(orbitalElements);

finalnu = table(nuFDeg,'VariableNames',"Final True Anomaly (deg)");
disp(finalnu);

keplerSME = table(KSME,'VariableNames',"Specific Mechanical Energy for Kepler's Method");
disp(keplerSME);

APKTable = table(APK,'VariableNames',"Specific Angular Momentum for Kepler's Method");
disp(APKTable);

%% Functions

function [a,e,i,RA,w,nu,mu] = RV2COE(R,V)
    mu = 398600;

    %Magnitude
    r = norm(R);
    v = norm(V);
    
    % Determination of Classical Orbital Elements
    % 1) Semi Major Axis
        epsillon = (v*v) / 2 - mu / r;
        a = - mu / (2 * epsillon);

    % 2) Eccentricity
        RVDot = dot(R,V);
        ebar = (1 / mu) * (R*(v*v-mu/r) - RVDot*V);
        e = norm(ebar);
    
    % 3) True Anomaly
        hbar = cross(R,V);
        h = norm(hbar);
        i = acosd(hbar(3)/h);
   
   % 4) Right Ascention of the Ascending Node
        kbar = [0,0,1];
        nbar = cross(kbar,hbar);
        n = norm(nbar);
        RA = acosd(nbar(1)/n);

        if nbar(2) < 0
            RA = 360 - RA;
        end

   % 5) Argument of Perigee
        nedot = dot(nbar,ebar);
        w = acosd(nedot/(n*e));
        if ebar(3) < 0
            w = 360 - w;
        end

   % 6) True Anomaly

        eRdot = dot(ebar, R);
        nu = acosd(eRdot/(e*r));   
        
        if RVDot < 0
            nu = 360 - nu;
        end


end

function [R,V] = COE2RV(a,e,i,RA,w,nu,mu)
    if(i == 0 || e == 0)
        disp("Warning!!! This is a special case! i or e = 0!")
    end
    P = a*(1-e*e);
    r = P / (1 + e*cosd(nu));

    snu = sind(nu);
    cnu = cosd(nu);

    Rpqw = [r*cnu; r*snu; 0];
    mup = sqrt(mu / P);
    Vpqw = [mup*(-snu); mup*(e+cnu); 0];
    
    cRA = cosd(RA);
    sRA = sind(RA);
    cw = cosd(w);
    sw = sind(w);
    ci = cosd(i);
    si = sind(i);
 
    A = [cRA*cw - sRA*sw*ci, -cRA*sw-sRA*cw*ci, sRA*si;
         sRA*cw+cRA*sw*ci, -sRA*sw+cRA*cw*ci, -cRA*si;
         sw*si, cw*si, ci
    ];

    R = A * Rpqw;
    V = A * Vpqw;

end
