function [ normed_vols ] = generate_vols( N, num_angles_1, num_angles_2, basis)
%GENERATE_VOLS Generates rotated vols
%   I'm going to genericize this with ability to add arbitrary numbers of
%   moving parts, etc.  But for right now, uses our good ol' benchmark

    if nargin > 3
        project_basis = 1;
    else
        project_basis = 0;
    end

    max_angle_1 = pi/2;
    max_angle_2 = pi/2;
    
    path_to_source_file = fileparts(mfilename('fullpath'));
    
    fixed_map = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapAlphaOne.mrc');
    vol_fixed = load_mrc(fixed_map);
    vol_fixed = max(0, vol_fixed);

    moving_map1 = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapAlphaTwo.mrc');
    vol_moving1 = load_mrc(moving_map1);
    vol_moving1 = max(0, vol_moving1);

    moving_map2 = fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapBeta.mrc');
    vol_moving2 = load_mrc(moving_map2);
    vol_moving2 = 0.85*max(0, vol_moving2);

    vol_moving1_rot = vol_rotate_z(vol_moving1, linspace(0, max_angle_1, num_angles_1));
    vol_moving2_rot = vol_rotate_z(vol_moving2, linspace(0, max_angle_2, num_angles_2));

    vol_moving_rot = bsxfun(@plus, vol_moving1_rot, permute(vol_moving2_rot, [1 2 3 5 4]));
    vol_moving_rot = reshape(vol_moving_rot, [size(vol_fixed) size(vol_moving1_rot, 4)*size(vol_moving2_rot, 4)]);

    vol_fixed = vol_downsample(vol_fixed, N);
    vol_moving_rot = vol_downsample(vol_moving_rot, N);

    vols = bsxfun(@plus, vol_fixed, vol_moving_rot);
    
    if project_basis
        vols = basis_project(basis, vols);
    end
    
%   we normalize the energy!
    avg_vol_energy = sum(vols(:).^2);
    if(avg_vol_energy ~= 0)
        normed_vols = vols/sqrt(avg_vol_energy) * sqrt(size(vols,4));
    else
        error('Vols have zero energy!');
    end
end
