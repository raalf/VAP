Fs = 1/valDELTIME;
T = 1/Fs;

L = size(Tip_Twist_Data,1);

t = (0:L-1)*T;

Y = fft(Tip_Twist_Data(:,end));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1)

% figure(10)
% hold on
% plot(valDELTIME*(valGUSTSTART:valTIMESTEP)-valGUSTSTART.*valDELTIME,(180/pi)*matTWISTGLOB(valGUSTSTART:end,end))
% ylabel('Tip Twist (deg)')
% grid on
% box on