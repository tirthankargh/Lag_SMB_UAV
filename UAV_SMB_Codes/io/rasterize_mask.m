function mask_raster = rasterize_mask(shapefile, R, X, Y)
    S = shaperead(shapefile);
    mask_raster = false(size(X));
    for k = 1:length(S)
        mask_raster = mask_raster | inpolygon(X, Y, S(k).X, S(k).Y);
    end
end
