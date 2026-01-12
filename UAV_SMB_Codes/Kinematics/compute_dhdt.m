function dhdt = compute_dhdt(DEM2, DEM1_adv, dt)
    dhdt = (DEM2 - DEM1_adv) / dt;
end
