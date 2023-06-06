function el = elevation(ENU, el_mask)

el=asin(ENU(:,3)./sqrt(ENU(:,1).^2+ENU(:,2).^2+ENU(:,3).^2)); %%Eelevation angle
el=el*180/pi; %%Radian -> Degree

%%계산된 elevation angle이 el_mask보다 커야만 정상적으로 출력하고, 크지 않을 경우 NaN으로 처리%%
a = size(el);
for i=1:a(1)
    if el(i)<el_mask
        el(i)= NaN;
    end
    
end

el=el';



