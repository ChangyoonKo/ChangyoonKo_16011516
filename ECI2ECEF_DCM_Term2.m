function DCM = ECI2ECEF_DCM_Term2(time,t)

w_earth = 360/86164; %[deg/sec] 

Theta_g = (siderealTime(juliandate(datetime(time)))+(w_earth*t*60))*pi/180; %% datetime->juliandate->siderealTime (Deg -> Radian)

DCM = [cos(Theta_g) sin(Theta_g) 0; -sin(Theta_g) cos(Theta_g) 0; 0 0 1;];