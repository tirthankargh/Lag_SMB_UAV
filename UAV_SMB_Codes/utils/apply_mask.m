function data_masked = apply_mask(data, mask)
    data_masked = data;
    data_masked(~mask) = NaN;
end
