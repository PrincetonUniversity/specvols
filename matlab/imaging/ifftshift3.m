function vol = ifftshift3(shifted_vol)
    vol = ifftshift(ifftshift(ifftshift(shifted_vol, 1), 2), 3);
end