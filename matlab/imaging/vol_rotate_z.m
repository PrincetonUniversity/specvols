function rotvol = vol_rotate_z(vol, angle_rad)
    rotvol = zeros(size(vol));
    nonzero_slices = find(squeeze(sum(sum(vol~=0, 1),2)));
    for i=1:length(nonzero_slices)
        z = nonzero_slices(i);
        rotvol(:,:,z) = imrotate(vol(:,:,z), radtodeg(angle_rad), 'bilinear', 'crop');
    end
end


