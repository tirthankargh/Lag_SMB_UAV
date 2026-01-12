function flux_div = compute_flux_divergence(H, vx, vy, f, res)
    % Compute ice flux in each direction
    % H: ice thickness (m)
    % vx, vy: horizontal velocity (m/yr)
    % f: vertical averaging factor (0-1), if not used set = 1
    % res: pixel resolution (m)
    

    % Apply factor
    ux = f .* vx;
    uy = f .* vy;

    % Ice flux
    Fx = H .* ux;  % m^2/yr
    Fy = H .* uy;  % m^2/yr

    % Compute flux divergence
    [dFx_dx, ~] = gradient(Fx, res); % ∂(Hx)/∂x
    [~, dFy_dy] = gradient(Fy, res); % ∂(Hy)/∂y

    flux_div = dFx_dx + dFy_dy;      % m/yr
end
