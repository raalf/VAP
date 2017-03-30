function [interp] = fcnPARABINTERP(vecY_KNOWN,vecX_KNOWN,vecX_DESIRED)

% This function performs a quadratic interpolation at the desired grid
% points, vecX_DESIRED, based on the known grid points, vecX_KNOWN, and
% known function values, vecY_KNOWN.
%
% NOTE: The vectors X_KNOWN and X_DESIRED MUST be in ascending order!

% Duplicate desired interpolation coordinates across columns to size of
% known coordinates
desired_points = repmat(vecX_DESIRED,1,size(vecX_KNOWN,1));

% Duplicate known coordinates across rows to size of interpolation
% coordinates
known_points = repmat(vecX_KNOWN',size(vecX_DESIRED,1),1);

% Difference 
diff = abs(known_points - desired_points);

% Determine the minimum differences in the matrix
min_loc = min(diff,[],2);

% Normalize the difference matrix by the minimum in each row
norm_min = diff./min_loc;

% Create index matrix for the location of the minimum values
idx = isnan(norm_min);
idx2 = (norm_min == 1);

idx = idx + idx2;

% Determine indices for a, b, and c coordinates
[~,idx_b] = find(idx == 1);

idx_a = idx_b - 1;
idx_c = idx_b + 1;

% Handle special case of indices at root of domain
idx_root = find(idx_a == 0);
idx_a(idx_root) = idx_a(idx_root) + 1;
idx_b(idx_root) = idx_b(idx_root) + 1;
idx_c(idx_root) = idx_c(idx_root) + 1;

% Handle special case of indices at tip of domain
idx_tip = find(idx_c == size(vecX_KNOWN,1)+1);
idx_a(idx_tip) = idx_a(idx_tip) - 1;
idx_b(idx_tip) = idx_b(idx_tip) - 1;
idx_c(idx_tip) = idx_c(idx_tip) - 1;

idx_x = (1:size(vecX_DESIRED,1))';

% Do quadratic interpolation
interp = vecY_KNOWN(idx_a)'.*(((vecX_DESIRED(idx_x) - vecX_KNOWN(idx_b)).*(vecX_DESIRED(idx_x) - vecX_KNOWN(idx_c)))./...
    ((vecX_KNOWN(idx_a)-vecX_KNOWN(idx_b)).*(vecX_KNOWN(idx_a) - vecX_KNOWN(idx_c)))) + vecY_KNOWN(idx_b)'.*...
    (((vecX_DESIRED(idx_x) - vecX_KNOWN(idx_c)).*(vecX_DESIRED(idx_x) - vecX_KNOWN(idx_a)))./...
    ((vecX_KNOWN(idx_b)-vecX_KNOWN(idx_c)).*(vecX_KNOWN(idx_b) - vecX_KNOWN(idx_a)))) +...
    vecY_KNOWN(idx_c)'.*(((vecX_DESIRED(idx_x) - vecX_KNOWN(idx_a)).*(vecX_DESIRED(idx_x) - vecX_KNOWN(idx_b)))./...
    ((vecX_KNOWN(idx_c)-vecX_KNOWN(idx_a)).*(vecX_KNOWN(idx_c) - vecX_KNOWN(idx_b))));

end