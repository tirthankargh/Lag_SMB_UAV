%% MAIN SMB PIPELINE
clear; clc; close all;
% Add subfolders
addpath('io','grids','kinematics','physics','utils');
% Load config
config;
% Make output folder
if ~exist(out_folder,'dir'); mkdir(out_folder); end

%% 1. Load input data
fprintf('Loading input data...\n');
% DEMs
[DEM1, R_DEM1] = load_geotiff(DEM1_file);
[DEM2, R_DEM2] = load_geotiff(DEM2_file);
DEM1 = double(DEM1);
DEM1(DEM1 < -1000) = NaN;
DEM2 = double(DEM2);
DEM2(DEM2 < -1000) = NaN;

info = geotiffinfo(DEM1_file);
proj = info.GeoTIFFTags.GeoKeyDirectoryTag;

% Ice thickness
[H, R_H] = load_geotiff(thickness_file);
H = double(H);
H(H < 0) = NaN;

% Velocity
[vx, R_vx] = load_geotiff(vx_file);
[vy, R_vy] = load_geotiff(vy_file);
vx = double(vx);
vx(abs(vx) > 1e4) = NaN;
vy = double(vy);
vy(abs(vy) > 1e4) = NaN;
%% 2. Define master 10 m grid from DEMs only
target_res = 10; % meters

% Use DEM extent only (not regional thickness/velocity products)
xmin = min([R_DEM1.XWorldLimits(1), R_DEM2.XWorldLimits(1)]);
xmax = max([R_DEM1.XWorldLimits(2), R_DEM2.XWorldLimits(2)]);
ymin = min([R_DEM1.YWorldLimits(1), R_DEM2.YWorldLimits(1)]);
ymax = max([R_DEM1.YWorldLimits(2), R_DEM2.YWorldLimits(2)]);

% Optional: add small buffer
buffer = 100; % 100 m buffer
xmin = xmin - buffer;
xmax = xmax + buffer;
ymin = ymin - buffer;
ymax = ymax + buffer;

fprintf('Glacier extent: %.2f x %.2f km\n', (xmax-xmin)/1000, (ymax-ymin)/1000);

xq = xmin : target_res : xmax;
yq = ymax : -target_res : ymin;  % descending to match MATLAB image coords

fprintf('Master grid: %d x %d pixels at %d m resolution\n', ...
        length(yq), length(xq), target_res);
fprintf('Estimated memory per array: %.1f MB\n', length(yq)*length(xq)*8/1e6);

% Safety check
if length(xq) * length(yq) > 100e6
    error('Grid too large (>100M pixels)! Check coordinate systems or increase target_res.');
end

[Xq, Yq] = meshgrid(xq, yq);

% Create spatial referencing object for master grid
R_master = maprefcells();

R_master.XWorldLimits = [min(xq) max(xq)];
R_master.YWorldLimits = [min(yq) max(yq)];
R_master.RasterSize   = size(Xq);
R_master.CellExtentInWorldX = target_res;
R_master.CellExtentInWorldY = target_res;
R_master.ColumnsStartFrom = 'north';
%% 3. Resample all datasets to master grid
fprintf('Resampling datasets to master grid...\n');
DEM1r = safe_resample(DEM1, R_DEM1, Xq, Yq);
DEM2r = safe_resample(DEM2, R_DEM2, Xq, Yq);
Hr    = safe_resample(H,   R_H,   Xq, Yq);
vxr   = safe_resample(vx,  R_vx,  Xq, Yq);
vyr   = safe_resample(vy,  R_vy,  Xq, Yq);
fprintf('Resampling complete.\n');

%% 4. Rasterize mask and apply
fprintf('Applying glacier mask...\n');
mask_r = rasterize_mask(mask_shapefile, R_master, Xq, Yq);  % returns logical array
DEM1r(~mask_r) = NaN;
DEM2r(~mask_r) = NaN;
Hr(~mask_r) = NaN;
vxr(~mask_r) = NaN;
vyr(~mask_r) = NaN;

%% 5. Flow-correct DEM1
fprintf('Flow-correcting DEM1...\n');
[DEM1_adv, adv_term] = flow_correct_dem(DEM1r, vxr, vyr, dt, Xq, Yq, target_res);

%% 6. Compute dh/dt
fprintf('Computing dh/dt...\n');
dhdt = compute_dhdt(DEM2r, DEM1_adv, dt);

%% 7. Compute flux divergence
fprintf('Computing flux divergence...\n');
f= 0.8 %column average velocity factor 
flux_div = compute_flux_divergence(Hr, vxr, vyr,f, target_res);

%% 8. Smooth flux divergence (choose method)
smoothing_option = 3; % 0=raw, 1=mean box, 2=gaussian, 3=exponential
fprintf('Smoothing flux divergence (method %d)...\n', smoothing_option);

switch smoothing_option
    case 0
        flux_div_smo = flux_div;
    case 1
        window_size = 3; % choose appropriate size
        flux_div_smo = smooth_flux_mean(flux_div, Hr, window_size);
    case 2
        sigma = 2; % standard deviation in pixels
        flux_div_smo = gauss_smooth(flux_div, sigma, Hr);
    case 3
        sigma = 2;
        alpha = 0.8;
        flux_div_smo = exp_smooth(flux_div, sigma, alpha, Hr);
end

%% 9. Compute SMB in m w.e./yr
fprintf('Computing SMB...\n');
SMB_ice = dhdt - adv_term + flux_div;
SMB = SMB_ice * rho_ice / rho_water;

%% DEBUG: Check data validity
fprintf('\n=== DATA VALIDITY CHECK ===\n');
fprintf('dhdt:        %d valid pixels (%.1f%%), range: [%.2f, %.2f] m/yr\n', ...
        sum(~isnan(dhdt(:))), 100*sum(~isnan(dhdt(:)))/numel(dhdt), ...
        min(dhdt(:)), max(dhdt(:)));
fprintf('flux_div_smo: %d valid pixels (%.1f%%), range: [%.2f, %.2f] m/yr\n', ...
        sum(~isnan(flux_div_smo(:))), 100*sum(~isnan(flux_div_smo(:)))/numel(flux_div_smo), ...
        min(flux_div_smo(:)), max(flux_div_smo(:)));
fprintf('SMB:         %d valid pixels (%.1f%%), range: [%.2f, %.2f] m w.e./yr\n', ...
        sum(~isnan(SMB(:))), 100*sum(~isnan(SMB(:)))/numel(SMB), ...
        min(SMB(:)), max(SMB(:)));
fprintf('Mask:        %d valid pixels (%.1f%%)\n', ...
        sum(mask_r(:)), 100*sum(mask_r(:))/numel(mask_r));

%% 10. Diagnostics
fprintf('\nGenerating diagnostic plots...\n');
%diagnostic_plots(dhdt, flux_div_smo, SMB);

figure('Position',[100 100 1500 500]);  % wide figure

% 1. dh/dt
subplot(1,3,1)
imagesc(dhdt); 
axis image; 
colormap(bluewhitered);  % or 'jet'
colorbar
title('dh/dt [m/yr]', 'FontSize',14)

% 2. Flux divergence
subplot(1,3,2)
imagesc(flux_div_smo); 
axis image; 
colormap(bluewhitered)
colorbar
title('Flux divergence [m/yr]', 'FontSize',14)

% 3. SMB
subplot(1,3,3)
imagesc(SMB); 
axis image; 
colormap(bluewhitered)
colorbar
title('SMB [m w.e./yr]', 'FontSize',14)
%% 11. Export results
fprintf('Exporting results...\n');
write_geotiff(fullfile(out_folder,'dhdt.tif'), dhdt, R_master, proj);
write_geotiff(fullfile(out_folder,'flux_div.tif'), flux_div_smo, R_master, proj);
write_geotiff(fullfile(out_folder,'SMB.tif'), SMB, R_master, proj);

% Save intermediate results for debugging (optional)
save(fullfile(out_folder,'intermediate_results.mat'), ...
     'DEM1_adv', 'flux_div', 'mask_r', 'Xq', 'Yq', 'R_master');

% Save metadata
metadata = struct();
metadata.target_resolution_m = target_res;
metadata.extent_km = [(xmax-xmin)/1000, (ymax-ymin)/1000];
metadata.grid_size = [length(yq), length(xq)];
metadata.dt_years = dt;
metadata.smoothing_method = smoothing_option;
metadata.processing_date = datestr(now);
save(fullfile(out_folder,'metadata.mat'), 'metadata');

fprintf('\n=== SMB pipeline completed successfully! ===\n');
fprintf('Output grid: %d x %d pixels (%.2f x %.2f km)\n', ...
        length(yq), length(xq), (xmax-xmin)/1000, (ymax-ymin)/1000);
fprintf('Results saved to: %s\n', out_folder);