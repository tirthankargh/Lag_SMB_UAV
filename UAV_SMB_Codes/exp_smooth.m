function flux_smo = exp_smooth(flux, sigma, alpha, H)
% Exponential smoothing of flux (approximated with Gaussian)
    % Apply Gaussian as a proxy for exponential kernel
    flux_smo = imgaussfilt(flux, sigma, 'FilterSize', 2*ceil(3*sigma)+1);
    % Apply alpha scaling
    flux_smo = alpha * flux_smo + (1-alpha) * flux;
    % Mass conservation
    flux_smo = flux_smo / (nansum(nansum(flux_smo)) / nansum(nansum(flux)));
end
