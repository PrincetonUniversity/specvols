% GEN_FAKEKV_VOLUMES generate the 3D volumes of the synthetic FakeKv model.
%
% This model consists of a central shaft with two blobs in at its ends that can rotate freely.
%
% Input:
%   L:
%       Side-length of the generated volumes.
%
%   max_angle_1, num_angles_1:
%       Take num_angles_1 rotations of the first blob, starting from zero with
%       angle delta equal to max_angle_1/num_angles_1.
%
%   max_angle_2, num_angles_2:
%       Same for the second blob.
%
% Output:
%   normed_vols
%       Normalized volumes, projected to Fourier-Bessel basis

function normed_vols = gen_fakekv_volumes(L, max_angle_1, num_angles_1, max_angle_2, num_angles_2)

    uncomputed_basis = fb_basis(L*ones(1,3),[]);
    basis = precompute_basis(uncomputed_basis);
    
    path_to_source_file = fileparts(mfilename('fullpath'));
    
    fixed_map = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapAlphaOne.mrc');
    vol_fixed = load_mrc(fixed_map);
    vol_fixed = max(0, vol_fixed);

    moving_map1 = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapAlphaTwo.mrc');
    vol_moving1 = load_mrc(moving_map1);
    vol_moving1 = max(0, vol_moving1);

    moving_map2 = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapBeta.mrc');
    vol_moving2 = load_mrc(moving_map2);
    vol_moving2 = 0.85*max(0, vol_moving2); % Make vol 2 darker -- Why?

    vol1rotations_downsampled = downsampled_z_rotations(L, vol_moving1, max_angle_1, num_angles_1);
    vol2rotations_downsampled = downsampled_z_rotations(L, vol_moving2, max_angle_2, num_angles_2);
    vol_fixed_downsampled = vol_downsample(vol_fixed, L);

    vols_preproj = zeros(L, L, L, num_angles_1, num_angles_2);
    for i=1:num_angles_1
        for j=1:num_angles_2
            vols_preproj(:,:,:,i,j) = vol1rotations_downsampled(:,:,:,i) + vol_fixed_downsampled + vol2rotations_downsampled(:,:,:,j);
        end
    end
    vols = basis_project(basis, reshape(vols_preproj, L,L,L,num_angles_1*num_angles_2));
    
    avg_vol_energy = sum(vols(:).^2);
    if(avg_vol_energy ~= 0)
        normed_vols = vols/sqrt(avg_vol_energy) * sqrt(size(vols,4));
    else
        error('Vols have zero energy!');
    end
end

function results = downsampled_z_rotations(L, volume, max_angle, num_angles)
    results = zeros(L, L, L, num_angles);
    for i=1:num_angles
        angle = -(i-1)*(max_angle/num_angles);
        rotated = vol_rotate_z(volume, angle);
        results(:,:,:,i) = vol_downsample(rotated, L);
    end
end
