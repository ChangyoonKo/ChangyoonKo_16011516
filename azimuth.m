function  az = azimuth(ENU)


az = atan2(ENU(:,1),ENU(:,2)); %%azimuth angle

%%부호 및 범위확인%%
a = size(az);
for i=1:a(1)
if az(i)<0
    az(i) = az(i)+2*pi;
end
end 

az = az'*180/pi; %%Radian -> Degree
