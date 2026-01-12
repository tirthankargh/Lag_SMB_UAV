function write_geotiff(filename, data, R, proj)
    geotiffwrite(filename, data, R, 'GeoKeyDirectoryTag', proj);
end
