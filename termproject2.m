clear; clc; close;

load nav.mat
GPS_a = nav.GPS.a*10^-3;          % semi-magor axis [km]
GPS_e = nav.GPS.e;                % Eccentrictiy
GPS_i = nav.GPS.i;                % Inclination [rad]
GPS_omega = nav.GPS.omega;        % Argument of perigee [rad]
GPS_Mean = nav.GPS.M0;            % Mean anomaly at toc M0 [rad]
GPS_toc = nav.GPS.toc;            % Epoch
GPS_OMEGA = nav.GPS.OMEGA;        % RAAN [rad]

QZSS_a = nav.QZSS.a*10^-3;        % semi-magor axis [km]
QZSS_e = nav.QZSS.e;              % Eccentrictiy
QZSS_i = nav.QZSS.i;              % Inclination [rad]
QZSS_omega = nav.QZSS.omega;      % Argument of perigee [rad]
QZSS_Mean = nav.QZSS.M0;          % Mean anomaly at toc M0 [rad]
QZSS_toc = nav.QZSS.toc;          % Epoch
QZSS_OMEGA = nav.QZSS.OMEGA;      % RAAN [rad]

BDS_a = nav.BDS.a*10^-3;          % semi-magor axis [km]
BDS_e = nav.BDS.e;                % Eccentrictiy
BDS_i = nav.BDS.i;                % Inclination [rad]
BDS_omega = nav.BDS.omega;        % Argument of perigee [rad]
BDS_Mean = nav.BDS.M0;            % Mean anomaly at toc M0 [rad]
BDS_toc = nav.BDS.toc;            % Epoch
BDS_OMEGA = nav.BDS.OMEGA;        % RAAN [rad]

mu = 3.986004418*10^5;            % [km^3 s^-2]
r2d = 180/pi;                     % Radian->Degree

GS_lat = 37;                      % Ground Station(Incheon) lat [deg]
GS_lon = 126;                     % Ground Station(Incheon) lon [deg]
GS_h = 1;                         % Ground Station h [km]
El_mask=0;                        % Elevation mask [deg]






for t=1:1:1440                    % 24h*60 minutes 1분마다 계산

    % GPS
    % PQW
    GPS_nu_PQW(t) = getnu(GPS_a,GPS_e,GPS_toc+[0 0 0 0 t 0],GPS_toc,GPS_Mean) * r2d;  %true_anomaly [degree]
    GPS_r_PQW(:,t) = solveRanglePerifocalFrame(GPS_a,GPS_e,GPS_nu_PQW(t));            % r of PQW
    GPS_v_PQW(:,t) = solveVelocityInPerifocalFrame(GPS_a,GPS_e,GPS_nu_PQW(t));        % v of PQW

    % PQW -> ECI
    GPS_r_ECI(:,t)= PQW2ECI(GPS_omega,GPS_i,GPS_OMEGA) * GPS_r_PQW(:,t);              % r of ECI
    GPS_v_ECI(:,t)= PQW2ECI(GPS_omega,GPS_i,GPS_OMEGA) * GPS_v_PQW(:,t);              % v of ECI

    % ECI -> ECEF
    GPS_r_ECEF(:,t) = ECI2ECEF_DCM_Term2(GPS_toc,t) * GPS_r_ECI(:,t);                 % r of ECEF
    GPS_v_ECEF(:,t) = ECI2ECEF_DCM_Term2(GPS_toc,t) * GPS_v_ECI(:,t);                 % v of ECEF

    % ECEF -> Geodetic
    wgs84 = wgs84Ellipsoid('kilometer');
    [a,b,c] = ecef2geodetic(wgs84,GPS_r_ECEF(1,t),GPS_r_ECEF(2,t),GPS_r_ECEF(3,t));
    GPS_lat(t)=a;
    GPS_lon(t)=b;
    GPS_h(t)=c;

    %ECEF -> ENU -> El,Az
    [GPS_E(t),GPS_N(t),GPS_U(t)] = ecef2enu(GPS_r_ECEF(1,t),GPS_r_ECEF(2,t),GPS_r_ECEF(3,t),GS_lat,GS_lon,GS_h,wgs84);
    GPS_ENU = [GPS_E(t),GPS_N(t),GPS_U(t)];
    GPS_El(t) = elevation(GPS_ENU,El_mask);                                            % [deg]
    GPS_Az(t) = azimuth(GPS_ENU);                                                      % [deg]




    % QZSS
    % PQW
    QZSS_nu_PQW(t) = getnu(QZSS_a,QZSS_e,QZSS_toc+[0 0 0 0 t 0],QZSS_toc,QZSS_Mean) * r2d;  %true_anomaly [degree]
    QZSS_r_PQW(:,t) = solveRanglePerifocalFrame(QZSS_a,QZSS_e,QZSS_nu_PQW(t));              % r of PQW
    QZSS_v_PQW(:,t) = solveVelocityInPerifocalFrame(QZSS_a,QZSS_e,QZSS_nu_PQW(t));          % v of PQW

    % PQW -> ECI
    QZSS_r_ECI(:,t)= PQW2ECI(QZSS_omega,QZSS_i,QZSS_OMEGA) * QZSS_r_PQW(:,t);               % r of ECI
    QZSS_v_ECI(:,t)= PQW2ECI(QZSS_omega,QZSS_i,QZSS_OMEGA) * QZSS_v_PQW(:,t);               % v of ECI

    % ECI -> ECEF
    QZSS_r_ECEF(:,t) = ECI2ECEF_DCM_Term2(QZSS_toc,t) * QZSS_r_ECI(:,t);                    % r of ECEF
    QZSS_v_ECEF(:,t) = ECI2ECEF_DCM_Term2(QZSS_toc,t) * QZSS_v_ECI(:,t);                    % v of ECEF

    % ECEF -> Geodetic
    wgs84 = wgs84Ellipsoid('kilometer');
    [a,b,c] = ecef2geodetic(wgs84,QZSS_r_ECEF(1,t),QZSS_r_ECEF(2,t),QZSS_r_ECEF(3,t));
    QZSS_lat(t)=a;
    QZSS_lon(t)=b;
    QZSS_h(t)=c;

    %ECEF -> ENU -> El,Az
    [QZSS_E(t),QZSS_N(t),QZSS_U(t)] = ecef2enu(QZSS_r_ECEF(1,t),QZSS_r_ECEF(2,t),QZSS_r_ECEF(3,t),GS_lat,GS_lon,GS_h,wgs84);
    QZSS_ENU = [QZSS_E(t),QZSS_N(t),QZSS_U(t)];
    QZSS_El(t) = elevation(QZSS_ENU,El_mask);                                             % [deg]
    QZSS_Az(t) = azimuth(QZSS_ENU);                                                       % [deg]

    % BDS
    % PQW
    BDS_nu_PQW(t) = getnu(BDS_a,BDS_e,BDS_toc+[0 0 0 0 t 0],BDS_toc,BDS_Mean) * r2d;      %true_anomaly [degree]
    BDS_r_PQW(:,t) = solveRanglePerifocalFrame(BDS_a,BDS_e,BDS_nu_PQW(t));                % r of PQW
    BDS_v_PQW(:,t) = solveVelocityInPerifocalFrame(BDS_a,BDS_e,BDS_nu_PQW(t));            % v of PQW

    % PQW -> ECI
    BDS_r_ECI(:,t)= PQW2ECI(BDS_omega,BDS_i,BDS_OMEGA) * BDS_r_PQW(:,t);                  % r of ECI
    BDS_v_ECI(:,t)= PQW2ECI(BDS_omega,BDS_i,GPS_OMEGA) * BDS_v_PQW(:,t);                  % v of ECI

    % ECI -> ECEF
    BDS_r_ECEF(:,t) = ECI2ECEF_DCM_Term2(BDS_toc,t) * BDS_r_ECI(:,t);                     % r of ECEF
    BDS_v_ECEF(:,t) = ECI2ECEF_DCM_Term2(BDS_toc,t) * BDS_v_ECI(:,t);                     % v of ECEF

    % ECEF -> Geodetic
    wgs84 = wgs84Ellipsoid('kilometer');
    [a,b,c] = ecef2geodetic(wgs84,BDS_r_ECEF(1,t),BDS_r_ECEF(2,t),BDS_r_ECEF(3,t));
    BDS_lat(t)=a;
    BDS_lon(t)=b;
    BDS_h(t)=c;

    %ECEF -> ENU -> El,Az
    [BDS_E(t),BDS_N(t),BDS_U(t)] = ecef2enu(BDS_r_ECEF(1,t),BDS_r_ECEF(2,t),BDS_r_ECEF(3,t),GS_lat,GS_lon,GS_h,wgs84);
    BDS_ENU = [BDS_E(t),BDS_N(t),BDS_U(t)];
    BDS_El(t) = elevation(BDS_ENU,El_mask);                                                % [deg]
    BDS_Az(t) = azimuth(BDS_ENU);                                                          % [deg]
    
end

%Plot Ground Track
figure(1)
geoplot(GPS_lat,GPS_lon,'r.')
hold on
geoplot(QZSS_lat,QZSS_lon,'g.')
hold on
geoplot(BDS_lat,BDS_lon,'b.')
legend('GPS','QZSS','BDS')

%Plot Sky View
figure(2)
skyplot(GPS_Az,GPS_El)
legend('GPS')
figure(3)
skyplot(QZSS_Az,QZSS_El)
legend('QZSS')
figure(4)
skyplot(BDS_Az,BDS_El)
legend('BDS')

%Plot satellite orbit
GPS_sc = satelliteScenario(datetime(GPS_toc),datetime(GPS_toc+[0 0 1 0 0 0]),60);
GPS_sat = satellite(GPS_sc,GPS_a*10^3,GPS_e,GPS_i*r2d,GPS_OMEGA*r2d,GPS_omega*r2d,GPS_nu_PQW(1));
GPS_GroundStation = groundStation(GPS_sc,GS_lat,GS_lon);
GPS_GroundTrack = groundTrack(GPS_sat);
satelliteScenarioViewer(GPS_sc)

QZSS_sc = satelliteScenario(datetime(QZSS_toc),datetime(QZSS_toc+[0 0 1 0 0 0]),60);
QZSS_sat = satellite(QZSS_sc,QZSS_a*10^3,QZSS_e,QZSS_i*r2d,QZSS_OMEGA*r2d,QZSS_omega*r2d,QZSS_nu_PQW(1));
QZSS_GroundStation = groundStation(QZSS_sc,GS_lat,GS_lon);
QZSS_GroundTrack = groundTrack(QZSS_sat);
satelliteScenarioViewer(QZSS_sc)

BDS_sc = satelliteScenario(datetime(BDS_toc),datetime(BDS_toc+[0 0 1 0 0 0]),60);
BDS_sat = satellite(BDS_sc,BDS_a*10^3,BDS_e,BDS_i*r2d,BDS_OMEGA*r2d,BDS_omega*r2d,BDS_nu_PQW(1));
BDS_GroundStation = groundStation(BDS_sc,GS_lat,GS_lon);
BDS_GroundTrack = groundTrack(BDS_sat);
satelliteScenarioViewer(BDS_sc)


