%% CONFIGURATION FILE
% Set file paths and physical constants

% File paths
DEM1_file      = '/Users/tirthankar/Research/Panchinala UAV/Inputs/DEM_PN_A_2021_clip_utm.tif';
DEM2_file      = '/Users/tirthankar/Research/Panchinala UAV/Inputs/PanchiA_Corrected_new_dsm.tif';
vx_file        = '/Users/tirthankar/Research/Panchinala UAV/Inputs/ITS_LIVE_reproject.tif';
vy_file        = '/Users/tirthankar/Research/Panchinala UAV/Inputs/ITS_LIVE_reproject_vy.tif';
thickness_file = '/Users/tirthankar/Research/Panchinala UAV/Inputs/THICKNESS_RGI-13-15.12_2022September22.tif';
mask_shapefile = '/Users/tirthankar/Research/Panchinala UAV/Panchinala/Panchi_mask.shp';

% Output folder
out_folder = '/Users/tirthankar/Research/Panchinala UAV/SMB_outputs';

% Physical constants
rho_ice = 900;      % kg/m3
rho_water = 1000;   % kg/m3

% Time difference between DEMs (years)
dt = 1;

% Vertical velocity scaling (vertically averaged factor)
F = 0.9;

% Target SMB resolution
target_res = 10; % meters

%% REPROJECT SHAPEFILE FROM LAT/LON TO UTM
fprintf('Reprojecting shapefile from geographic to UTM coordinates...\n');

% Read original shapefile (in lat/lon WGS84)
S_latlon = shaperead(mask_shapefile);

% Define projection systems
% WGS84 Geographic (EPSG:4326)
geo_crs = geocrs(4326);

% UTM Zone 43N (EPSG:32643) - for 77°E longitude in Northern Hemisphere
utm_crs = projcrs(32643);

% Reproject each feature in the shapefile
S_utm = S_latlon;
for k = 1:length(S_utm)
    % Get lat/lon coordinates
    lat = S_latlon(k).Y;
    lon = S_latlon(k).X;
    
    % Find NaN separators
    valid_idx = ~isnan(lat) & ~isnan(lon);
    
    % Reproject valid coordinates
    if any(valid_idx)
        [x_utm, y_utm] = projfwd(utm_crs, lat(valid_idx), lon(valid_idx));
        
        % Put back into structure with NaN separators in same positions
        S_utm(k).X = nan(size(lon));
        S_utm(k).Y = nan(size(lat));
        S_utm(k).X(valid_idx) = x_utm;
        S_utm(k).Y(valid_idx) = y_utm;
    end
    
    % Update geometry type if needed
    if isfield(S_utm(k), 'Geometry')
        S_utm(k).Geometry = S_latlon(k).Geometry;
    end
end

% Save reprojected shapefile (optional - for future use)
mask_shapefile_utm = fullfile(out_folder, 'Panchi_mask_utm.shp');
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end

try
    shapewrite(S_utm, mask_shapefile_utm);
    fprintf('Reprojected shapefile saved to: %s\n', mask_shapefile_utm);
    % Use the reprojected shapefile
    mask_shapefile = mask_shapefile_utm;
catch
    warning('Could not save reprojected shapefile, will use in-memory version');
    % Store in a variable that can be used instead
    mask_shapefile_data = S_utm;
end

% Verify projection
fprintf('Original shapefile extent: X=[%.2f, %.2f]°, Y=[%.2f, %.2f]°\n', ...
        min([S_latlon.X]), max([S_latlon.X]), min([S_latlon.Y]), max([S_latlon.Y]));
fprintf('Reprojected extent: X=[%.2f, %.2f] m, Y=[%.2f, %.2f] m\n', ...
        min([S_utm.X]), max([S_utm.X]), min([S_utm.Y]), max([S_utm.Y]));

fprintf('Configuration loaded successfully!\n\n');