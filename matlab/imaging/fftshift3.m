function shifted_vol = fftshift3(vol)
    shifted_vol = fftshift(fftshift(fftshift(vol, 1), 2), 3);
end