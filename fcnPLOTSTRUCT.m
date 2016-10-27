function [hFig4, hFig5] = fcnPLOTSTRUCT(matDEFLECTION, matTWIST)

% Get total number of timesteps stored in matDEFLECTION
nTIMESTEPS = rows(matDEFLECTION);

hFig4 = figure(4);

% Plot deflection at latest timestep
plot(matDEFLECTION(nTIMESTEPS,:))

xlabel('Spanwise location (m)');
ylabel('Wing Deflection (m)');
title('Dynamic Loading - Deflection');

hFig5 = figure(5);

% Plot twist at latest timestep
plot(matTWIST(nTIMESTEPS,:))

xlabel('Spanwise location (m)');
ylabel('Wing Twist (deg)');
title('Dynamic Loading - Twist');

end
