function rangelnPQW = solveRanglePerifocalFrame(semimajor_axis,eccentricity,true_anomaly) %%'km' 'deg'


p = semimajor_axis*(1-(eccentricity)^2);
nu = true_anomaly*(pi/180); %%degree to radian
r = p/(1+eccentricity*cos(nu));
rangelnPQW = [r*cos(nu); r*sin(nu); 0];