function DCM = ECI2ECEF_DCM(time)

Theta_g = siderealTime(juliandate(datetime(time)))*pi/180; %% datetime->juliandate->siderealTime (Deg -> Radian)

DCM = [cos(Theta_g) sin(Theta_g) 0; -sin(Theta_g) cos(Theta_g) 0; 0 0 1;];