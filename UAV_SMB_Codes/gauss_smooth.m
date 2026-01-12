function flux_smo = gauss_smooth(flux, sigma, H)
% Gaussian smoothing of flux
    flux_smo = imgaussfilt(flux, sigma, 'FilterSize', 2*ceil(3*sigma)+1);
    % mass-conservation normalization
    flux_smo = flux_smo / (nansum(nansum(flux_smo)) / nansum(nansum(flux)));
end
