function [DEM1_adv, adv_term] = flow_correct_dem(DEM1, vx, vy, dt, Xq, Yq, res)
    % Forward/backward project DEM1 using velocities
    % Also computes slope-parallel advection term: us · ∇h

    % Compute displacement
    dx = vx * dt;
    dy = vy * dt;

    % Backward advection to DEM2 coordinates
    Xb = Xq - dx;
    Yb = Yq - dy;

    % Interpolate DEM1 to new positions
    DEM1_adv = interp2(Xq, Yq, DEM1, Xb, Yb, 'linear', NaN);

    % Compute slope (gradients) of original DEM1
    [dh_dx, dh_dy] = gradient(DEM1, res);

    % Slope-parallel advection
    adv_term = vx .* dh_dx + vy .* dh_dy;  % m/yr
end
