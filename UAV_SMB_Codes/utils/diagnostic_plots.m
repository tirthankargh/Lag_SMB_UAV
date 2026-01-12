function diagnostic_plots(dhdt, flux_div, SMB)
    figure; imagesc(dhdt); colormap(bluewhitered); colorbar; title('dh/dt (m/yr)'); axis equal tight
    figure; imagesc(flux_div); colormap(bluewhitered); colorbar; title('Flux divergence (m/yr)'); axis equal tight
    figure; imagesc(SMB); colormap(bluewhitered); colorbar; title('SMB (m w.e./yr)'); axis equal tight
end
