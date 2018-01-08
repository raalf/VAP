valGUSTL = 6;
valUINF = 50;
valGUSTAMP = 4;
gust_vel = 0;
matUINF = [50 0 0];
i = 1;

for delx = 0:0.001:30
    
    idx3 = delx(:,1) >= 0 & delx(:,1) < 0.25*valGUSTL;
    idx3_1 = find(idx3 > 0);
    idx4 = delx(:,1) >= 0.25*valGUSTL & delx(:,1) < 0.75*valGUSTL;
    idx4_1 = find(idx4 > 0);
    idx5 = delx(:,1) >= 0.75*valGUSTL & delx(:,1) <= valGUSTL;
    idx5_1 = find(idx5 > 0);

    if any(idx3 > 0)
        tau = delx(idx3_1,1)./valUINF;
        gust_vel(idx3_1) = 0.5*valGUSTAMP*(1 - cos((2*pi*tau/(0.5*valGUSTL/valUINF))));
        matUINF(idx3_1,3) = matUINF(idx3_1,3) + (gust_vel(idx3_1));
    end

    if any(idx4 > 0)
        tau = delx(idx4_1,1)./valUINF;
        gust_vel(idx4_1) = -valGUSTAMP*(cos(((2*pi/(valGUSTL/valUINF)))*tau + (3*pi/(2*valGUSTL/valUINF))));
        matUINF(idx4_1,3) = matUINF(idx4_1,3) + (gust_vel(idx4_1));
    end

    if any(idx5 > 0)
        tau = delx(idx5_1,1)./valUINF;
        gust_vel(idx5_1) = -0.5*valGUSTAMP*(1 - cos((2*pi*tau/(0.5*valGUSTL/valUINF))));
        matUINF(idx5_1,3) = matUINF(idx5_1,3) + (gust_vel(idx5_1));
    end
    
    x(i) = gust_vel;
    i = i + 1;
end

figure(1)
plot((0:0.001:30)./valUINF,x)

% gust_vel = 0;
% matUINF = [50 0 0];
% i = 1;
% 
% for delx = 0:0.01:6
%     
%     idx3 = delx(:,1) >= 0 & delx(:,1) < 0.5*valGUSTL;
%     idx3_1 = find(idx3 > 0);
%     idx4 = delx(:,1) >= 0.5*valGUSTL & delx(:,1) <= valGUSTL;
%     idx4_1 = find(idx4 > 0);
% %     idx5 = delx(:,1) >= 0.75*valGUSTL & delx(:,1) <= valGUSTL;
% %     idx5_1 = find(idx5 > 0);
% 
%     if any(idx3 > 0)
%         tau = delx(idx3_1,1)./valUINF;
%         gust_vel(idx3_1) = 0.5*valGUSTAMP*(1 - cos((2*pi*tau/(0.5*valGUSTL/valUINF))));
%         matUINF(idx3_1,3) = matUINF(idx3_1,3) + (gust_vel(idx3_1));
%     end
% 
%     if any(idx4 > 0)
%         tau = delx(idx4_1,1)./valUINF;
%         gust_vel(idx4_1) = 0.5*valGUSTAMP*(-1+cos(((2*pi/(0.5*valGUSTL/valUINF)))*tau));
%         matUINF(idx4_1,3) = matUINF(idx4_1,3) + (gust_vel(idx4_1));
%     end
% 
% %     if any(idx5 > 0)
% %         tau = delx(idx5_1,1)./valUINF;
% %         gust_vel(idx5_1) = -0.5*valGUSTAMP*(1 - cos((2*pi*tau/(0.5*valGUSTL/valUINF))));
% %         matUINF(idx5_1,3) = matUINF(idx5_1,3) + (gust_vel(idx5_1));
% %     end
%     
%     y(i) = gust_vel;
%     i = i + 1;
% end
% 
% figure(2)
% plot(0:0.01:6,y)
