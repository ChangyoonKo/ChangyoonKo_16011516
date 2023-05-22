function PQW2ECI = PQW2ECI(arg_prg,inc_angle,RAAN)

a=[cos(arg_prg) -sin(arg_prg) 0; sin(arg_prg) cos(arg_prg) 0; 0 0 1]; %%R(-W,3")
b=[1 0 0; 0 cos(inc_angle) -sin(inc_angle); 0 sin(inc_angle) cos(inc_angle)];%%R(-i,1)
c=[cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];%%R(-Ω,3)
PQW2ECI=a*b*c;
