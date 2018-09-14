% Example parameters:
% N = 8
% n = 2^9 - 1
% max_angle_1 = pi/2;
% num_angles_1 = 16;
% max_angle_2 = 0;
% num_angles_2 = 1;

% GEN_FAKEKV_VOLUMES generate the 3D volumes of the synthetic FakeKv model.
%
% This model consists of a central shaft with two blobs in at its ends that can rotate freely.
%
% Input:
%   L:
%       Side-length of the generated volumes.
%
%   max_angle_1/num_angles_1:
%       Take num_angles_1 rotations of the first blob starting from zero
%       with angle delta equal to max_angle_1/num_angles_2
%
%   max_angle_2/num_angles_2:
%       Same for the second blob.
%
% Output:
%   normed_vols

function normed_vols = gen_fakekv_volumes(L, max_angle_1, num_angles_1, max_angle_2, num_angles_2)

    uncomputed_basis = fb_basis(L*ones(1,3),[]);
    basis = precompute_basis(uncomputed_basis);

    path_to_source_file = fileparts(mfilename('fullpath'));
    fixed_map = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapAlphaOne.mrc');
    moving_map1 = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapAlphaTwo.mrc');
    moving_map2 = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapBeta.mrc');
    vol_fixed = load_mrc(fixed_map);
    vol_moving1 = load_mrc(moving_map1);
    vol_moving2 = load_mrc(moving_map2);

    vol_fixed = max(0, vol_fixed);
    vol_moving1 = max(0, vol_moving1);
    vol_moving2 = 0.85*max(0, vol_moving2);

    angles_1 = linspace(0, max_angle_1, num_angles_1);
    angles_2 = linspace(0, max_angle_2, num_angles_2);

    rots_1 = angles_to_rots(cat(1, angles_1, zeros(2, num_angles_1)));
    rots_2 = angles_to_rots(cat(1, angles_2, zeros(2, num_angles_2)));
    
    vol_moving1_rot = vol_rotate(vol_moving1, rots_1);
    vol_moving2_rot = vol_rotate(vol_moving2, rots_2);

    vol_moving1_rot_downsampled = vol_downsample(vol_moving1_rot, L);
    vol_moving2_rot_downsampled = vol_downsample(vol_moving2_rot, L);
    
    vol_fixed_downsampled = vol_downsample(vol_fixed, L);

    vol_moving_both_rot_downsampled = bsxfun(@plus, vol_moving1_rot_downsampled, permute(vol_moving2_rot_downsampled, [1 2 3 5 4]));
    vol_moving_both_rot_downsampled = reshape(vol_moving_both_rot_downsampled, [size(vol_fixed_downsampled) size(vol_moving1_rot_downsampled, 4)*size(vol_moving2_rot_downsampled, 4)]);
    
    vols_preproj = bsxfun(@plus, vol_fixed_downsampled, vol_moving_both_rot_downsampled);
    vols = basis_project(basis, vols_preproj);

    %we'll normalize anyway!
    avg_vol_energy = sum(vols(:).^2);
    if(avg_vol_energy ~= 0)
        normed_vols = vols/sqrt(avg_vol_energy) * sqrt(size(vols,4));
    else
        error('Vols have zero energy!');
    end
end
