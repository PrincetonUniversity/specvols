% GEN_FAKEKV_TOP_ROTATIONS generate rotations of the top part of the
% synthetic FakeKV molecule
%
% This model consists of a central shaft with two blobs in at its ends that can rotate freely.
%
% Input:
%   L:
%       Side-length of the generated volumes.
%
%   angles:
%       Vector of angles (rad)
%
% Output:
%   rotated_top_vols
%       rotated top volumes, in the standard basis, not normalized

function rotated_top_vols = gen_fakekv_top_rotations(L, angles)
    path_to_source_file = fileparts(mfilename('fullpath'));

    vol_top = max(0, load_mrc(fullfile(path_to_source_file, '../../data/FakeKvMap/FakeKvMapAlphaTwo.mrc')));

    rotated_top_vols = zeros(L, L, L, length(angles));
    for i=1:length(angles)
        rotated = vol_rotate_z(vol_top, angles(i));
        rotated_top_vols(:,:,:,i) = vol_downsample(rotated, L);
    end
end

