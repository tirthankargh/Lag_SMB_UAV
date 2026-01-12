function [data, R] = load_geotiff(filename)
    [data, R] = readgeoraster(filename);
    data = double(data);
end
