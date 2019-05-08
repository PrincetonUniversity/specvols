function results = downsampled_z_rotations(L, volume, max_angle, num_angles)
    results = zeros(L, L, L, num_angles);
    for i=1:num_angles
        angle = -(i-1)*(max_angle/num_angles);
        rotated = vol_rotate_z(volume, angle);
        results(:,:,:,i) = vol_downsample(rotated, L);
    end
end
