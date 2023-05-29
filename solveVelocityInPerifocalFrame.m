function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis,eccentricity,true_anomaly) %%'km' 'deg'

mu = 3.986004418*10^5; %%[km^3 s^-2]
p = semimajor_axis*(1-(eccentricity)^2);
nu = true_anomaly*(pi/180); %%degree to radian
velocityInPQW = sqrt(mu/p)*[-sin(nu); eccentricity+cos(nu); 0];

