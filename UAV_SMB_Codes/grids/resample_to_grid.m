function data_resamp = resample_to_grid(data, R, Xq, Yq)
    % Ensure double and NaN for no-data
    data = double(data);
    data(data <= -9999) = NaN;  % handles all your no-data placeholders

    [nrows, ncols] = size(data);

    % Pixel center coordinates
    x_vec = R.XWorldLimits(1) + R.CellExtentInWorldX/2 : R.CellExtentInWorldX : R.XWorldLimits(2) - R.CellExtentInWorldX/2;
    y_vec = R.YWorldLimits(2) - R.CellExtentInWorldY/2 : -R.CellExtentInWorldY : R.YWorldLimits(1) + R.CellExtentInWorldY/2;

    % Safety check
    assert(length(x_vec) == ncols, 'X vector length mismatch!');
    assert(length(y_vec) == nrows, 'Y vector length mismatch!');

    % Make meshgrid
    [X, Y] = meshgrid(x_vec, y_vec);

    % Interpolate to master grid
    F = griddedInterpolant(X, Y, data, 'linear', NaN);
    data_resamp = F(Xq, Yq);
end
