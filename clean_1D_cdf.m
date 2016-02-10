function x = clean_1D_cdf(c)
    c.d1 = [];
    c.d2 = [];
    c.scan_index = [];
    c.tile_idx = [];
    x = c.xyz; 
    clearvars c
end

