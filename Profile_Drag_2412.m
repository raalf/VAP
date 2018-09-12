clear cd cd_t sectionCL
% load naca2412Mike.dat

valTIMEPER = ceil((valGUSTL/valUINF)/valDELTIME); % Number of time steps for one gust period

int_stop = 3*valTIMEPER + valGUSTSTART; % Time step to stop time integration

CD0 = vecCD(valGUSTSTART); % Reference total drag coefficient
CD = mean(vecCD(valGUSTSTART:int_stop));

integrand_t = (vecCD(valGUSTSTART:int_stop)-CD0)./CD0;

reduction = (CD - CD0)/CD0