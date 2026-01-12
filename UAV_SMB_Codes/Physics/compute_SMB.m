% 1. Compute flux divergence
flux_div = compute_flux_divergence(Hr, vxr, vyr, f, target_res);

% 2. Compute slope-parallel advection
adv_term = compute_advection(DEM1r, vxr, vyr, target_res);

% 3. SMB in ice equivalent (m/yr)
SMB_ice = dhdt - adv_term + flux_div;  

% 4. Convert to water equivalent
rho_ice = 900; % kg/m3
rho_water = 1000; % kg/m3
SMB = SMB_ice * rho_ice / rho_water; % m w.e./yr
