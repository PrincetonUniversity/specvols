function shifted_vol = fakekv_bottom_shift(vol, xshift, yshift)
    assert(isequal(size(vol), [108 108 108]));
    L = 108;
    
    [X,Y,Z] = meshgrid(1:L, 1:L, 1:L);
    shifted_vol = zeros(L,L,L);
    
    Z_START = 16;
    Z_END = 54;
    
    for z = Z_START:Z_END
        shift_amount = ((Z_END - z)/(Z_END-Z_START))^2;
        shifted_vol(:,:,z) = shift_image(vol(:,:,z), xshift*shift_amount, yshift*shift_amount, 0);
    end
    shifted_vol(:,:,1:Z_START-1) = vol(:,:,1:Z_START-1);
    shifted_vol(:,:,Z_END+1:end) = vol(:,:,Z_END+1:end);
    shifted_vol(isnan(shifted_vol)) = 0;
end