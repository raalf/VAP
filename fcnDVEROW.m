function [matROWS] = fcnDVEROW(ledves, vecDVEPANEL, vecDVEWING, vecM, vecN)

% This functions determines what DVEs are in which row or column across the
% span of the wing. matROWS will be a matrix that is NELE x vecM.

lepanels = vecDVEPANEL(ledves);

% Determine DVEs in each spanwise station
for i = 1:max(vecDVEWING)

	idxdve = ledves(vecDVEWING(ledves) == i);
	idxpanel = lepanels(vecDVEWING(ledves) == i);

    m = vecM(idxpanel);
    if any(m - m(1))
        disp('Problem with wing chordwise elements.');
        break
    end
    m = m(1);

    tempm = repmat(vecN(idxpanel), 1, m).*repmat([0:m-1],length(idxpanel),1);
    
    matROWS = repmat(idxdve,1,m) + tempm;

end

end