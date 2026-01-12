function data_resamp = safe_resample(data, R, Xq, Yq)
    % Convert to double and clean data
    data = double(data);
    data(~isfinite(data)) = NaN;
    
    [nrows, ncols] = size(data);
    
    fprintf('    Input data: %dx%d, %d valid pixels\n', nrows, ncols, sum(~isnan(data(:))));
    
    % Handle both MapCellsReference and MapPostingsReference
    if isa(R, 'map.rasterref.MapCellsReference')
        % Cells: pixel centers
        x_vec = R.XWorldLimits(1) + R.CellExtentInWorldX/2 : ...
                R.CellExtentInWorldX : ...
                R.XWorldLimits(2) - R.CellExtentInWorldX/2;
        y_vec = R.YWorldLimits(2) - R.CellExtentInWorldY/2 : ...
                -R.CellExtentInWorldY : ...
                R.YWorldLimits(1) + R.CellExtentInWorldY/2;
    else
        % Postings: use linspace for evenly-spaced coordinates
        x_vec = linspace(R.XWorldLimits(1), R.XWorldLimits(2), ncols);
        y_vec = linspace(R.YWorldLimits(2), R.YWorldLimits(1), nrows);
    end
    
    fprintf('    Coordinate vectors: x_vec has %d points, y_vec has %d points\n', ...
            length(x_vec), length(y_vec));
    fprintf('    x_vec range: [%.2f, %.2f]\n', min(x_vec), max(x_vec));
    fprintf('    y_vec range: [%.2f, %.2f]\n', min(y_vec), max(y_vec));
    
    % Resize raster if X/Y mismatch
    if length(x_vec) ~= ncols || length(y_vec) ~= nrows
        fprintf('    WARNING: Resizing data from %dx%d to %dx%d\n', ...
                nrows, ncols, length(y_vec), length(x_vec));
        data = imresize(data, [length(y_vec), length(x_vec)], 'bilinear');
    end
    
    % Remove rows and columns that are all NaN
    valid_rows = ~all(isnan(data), 2);
    valid_cols = ~all(isnan(data), 1);
    
    fprintf('    Valid rows: %d/%d, Valid cols: %d/%d\n', ...
            sum(valid_rows), length(valid_rows), sum(valid_cols), length(valid_cols));
    
    % Check if any valid data exists
    if ~any(valid_rows) || ~any(valid_cols)
        warning('No valid data found, returning NaN grid');
        data_resamp = NaN(size(Xq));
        return;
    end
    
    % Extract valid data and coordinate vectors
    data = data(valid_rows, valid_cols);
    x_vec = x_vec(valid_cols);
    y_vec = y_vec(valid_rows);
    
    fprintf('    After removing NaN rows/cols: data is %dx%d\n', size(data,1), size(data,2));
    fprintf('    After trimming: x_vec has %d points [%.2f, %.2f]\n', ...
            length(x_vec), min(x_vec), max(x_vec));
    fprintf('    After trimming: y_vec has %d points [%.2f, %.2f]\n', ...
            length(y_vec), min(y_vec), max(y_vec));
    
    % griddedInterpolant requires ASCENDING order
    if y_vec(1) > y_vec(end)
        fprintf('    Flipping y_vec to ascending order\n');
        y_vec = flip(y_vec);
        data = flipud(data);
    end
    
    if x_vec(1) > x_vec(end)
        fprintf('    Flipping x_vec to ascending order\n');
        x_vec = flip(x_vec);
        data = fliplr(data);
    end
    
    fprintf('    Final x_vec: [%.2f, %.2f], y_vec: [%.2f, %.2f]\n', ...
            min(x_vec), max(x_vec), min(y_vec), max(y_vec));
    fprintf('    Query Xq: [%.2f, %.2f], Yq: [%.2f, %.2f]\n', ...
            min(Xq(:)), max(Xq(:)), min(Yq(:)), max(Yq(:)));
    
    % Check overlap
    x_overlap = (max(Xq(:)) >= min(x_vec)) && (min(Xq(:)) <= max(x_vec));
    y_overlap = (max(Yq(:)) >= min(y_vec)) && (min(Yq(:)) <= max(y_vec));
    fprintf('    X overlap: %s, Y overlap: %s\n', mat2str(x_overlap), mat2str(y_overlap));
    
    if ~x_overlap || ~y_overlap
        warning('Query grid does not overlap with data grid!');
        data_resamp = NaN(size(Xq));
        return;
    end
    
    % Use vectors directly for griddedInterpolant
    try
        F = griddedInterpolant({y_vec, x_vec}, data, 'linear', 'none');
        data_resamp = F(Yq, Xq);
        fprintf('    Interpolation successful: %d valid output pixels\n', sum(~isnan(data_resamp(:))));
    catch ME
        warning('griddedInterpolant failed: %s', ME.message);
        data_resamp = NaN(size(Xq));
    end
end