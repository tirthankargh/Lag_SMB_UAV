%% MAIN SMB PIPELINE
clear; clc; close all;

% Add subfolders
addpath('io','grids','kinematics','physics','utils');

% Load config
config;

% Make output folder
if ~exist(out_folder,'dir'); mkdir(out_folder); end

%% 1. Load input data
% DEMs
[DEM1, R_DEM1] = load_geotiff(DEM1_file);
[DEM2, R_DEM2] = load_geotiff(DEM2_file);
DEM1 = double(DEM1);
DEM1(DEM1 < -1000) = NaN;
DEM2 = double(DEM2);
DEM2(DEM2 < -1000) = NaN;

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

%% 2. Define master 10 m grid
target_res = 10; % meters
xmin = min([R_DEM1.XWorldLimits(1), R_DEM2.XWorldLimits(1)]);
xmax = max([R_DEM1.XWorldLimits(2), R_DEM2.XWorldLimits(2)]);
ymin = min([R_DEM1.YWorldLimits(1), R_DEM2.YWorldLimits(1)]);
ymax = max([R_DEM1.YWorldLimits(2), R_DEM2.YWorldLimits(2)]);

xq = xmin : target_res : xmax;
yq = ymax : -target_res : ymin;  % descending to match MATLAB image coords
[Xq, Yq] = meshgrid(xq, yq);

%% 3. Resample all datasets to master grid
DEM1r = safe_resample(DEM1, R_DEM1, Xq, Yq);
DEM2r = safe_resample(DEM2, R_DEM2, Xq, Yq);
Hr    = safe_resample(H,   R_H,   Xq, Yq);
vxr   = safe_resample(vx,  R_vx,  Xq, Yq);
vyr   = safe_resample(vy,  R_vy,  Xq, Yq);

%% 4. Rasterize mask and apply
mask_r = rasterize_mask(mask_shapefile, Xq, Yq);  % returns logical array
DEM1r(~mask_r) = NaN; 
DEM2r(~mask_r) = NaN;
Hr(~mask_r) = NaN; 
vxr(~mask_r) = NaN; 
vyr(~mask_r) = NaN;

%% 5. Flow-correct DEM1
DEM1_adv = flow_correct_dem(DEM1r, vxr, vyr, dt, Xq, Yq, target_res);

%% 6. Compute dh/dt
dhdt = compute_dhdt(DEM2r, DEM1_adv, dt);

%% 7. Compute flux divergence
flux_div = compute_flux_divergence(Hr, vxr, vyr, F, target_res);

Smooth flux divergence (choose method)
smoothing_option = 2; % 0=raw, 1=mean box, 2=gaussian, 3=exponential

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
SMB = compute_SMB(dhdt, flux_div_smo, rho_ice, rho_water);

%% Diagnostics
diagnostic_plots(dhdt, flux_div, SMB);

%% Export results
write_geotiff(fullfile(out_folder,'dhdt.tif'), dhdt, R_DEM1);
write_geotiff(fullfile(out_folder,'flux_div.tif'), flux_div, R_DEM1);
write_geotiff(fullfile(out_folder,'SMB.tif'), SMB, R_DEM1);

disp('SMB pipeline completed!');
