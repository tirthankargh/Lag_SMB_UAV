function flux_smo = smooth_flux_mean(flux, H, window_size)
% Mean box smoothing of flux
    kernel = ones(window_size) / window_size^2;
    flux_smo = conv2(flux, kernel, 'same');
    % mass-conservation normalization
    flux_smo = flux_smo / (nansum(nansum(flux_smo)) / nansum(nansum(flux)));
end
