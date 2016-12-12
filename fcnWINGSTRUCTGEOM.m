function [vecSPNWSECRD, vecSPNWSEAREA, matQTRCRD, vecQTRCRD] = fcnWINGSTRUCTGEOM(vecDVEWING, vecDVELE, vecDVEPANEL, vecM, vecN, vecDVEHVCRD, matDVE, matVLST, vecDVEAREA)
    
    vecSPNWSECRD = [];
    vecSPNWSEAREA = [];
    [ledves, ~, ~] = find(vecDVELE > 0);
    lepanels = vecDVEPANEL(ledves);

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

        rows = repmat(idxdve,1,m) + tempm; % DVEs along each chord station

       vecSPNWSECRD = [vecSPNWSECRD; 2*sum(vecDVEHVCRD(rows),2)];  % Chord length at each spanwise station mid-point
       vecSPNWSEAREA = [vecSPNWSEAREA; sum(vecDVEAREA(rows),2)];

    end
    
    tempMIDLE = matVLST(matDVE(:,1:2));
    vecMIDLE = (tempMIDLE(:,1) + tempMIDLE(:,2))./2; % X location of mid-point on each DVE LE
    
    matMIDLE(:,1:vecM(1)) = vecMIDLE(rows); % DVE LE mid-point location at each chord station
    
    vecQTRCRD = matMIDLE(:,1) + vecSPNWSECRD*0.25; % Aerodynamic center line
    matQTRCRD = -(matMIDLE - repmat(vecQTRCRD,1,vecM(1))); % Distance from DVE LE mid-point to AC

end