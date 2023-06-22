function nu = getnu(a,e,t,toc,M0)

mu = 3.986004418*10^5; %%[km^3 s^-2]

dt = datetime(t)-datetime(toc);

[h,m,s] = hms(dt);

ds = h*3600+m*60+s;

n = sqrt(mu/(a^3));

M = n*ds + M0;

E = M;

while(1)

    E1 = E;

    E = M+e*sin(E);
    if abs(E-E1)<=0.00001
        break
    end
end

nu = 2*atan2(tan(E/2)*sqrt(1+e),sqrt(1-e));

if nu<0
    nu = nu+2*pi;
end


