% GEN_FAKEKV_MID_BOTTOM_SHIFTS generate shifts or strechings along the x-y
% plane of the middle and bottom part of the FakeKV molecule.
%
% Input:
%   L:
%       Side-length of the generated volumes.
%
%   xyshifts:
%       2 by n vector of shifts (in pixels) of the bottom-most slice of
%       of the bottom of the FakeKV molecule. The rest of the X-Y slices,
%       starting from the middle of the volume are also shifted in
%       a continuous manner.
%
% Output:
%   shifted_vols
%       middle and bottom parts of the FakeKV, stretched,
%       in the standard basis, not normalized.

function shifted_vols = gen_fakekv_mid_bottom_shifts(L, xyshifts)
    [two,n] = size(xyshifts);
    assert(two == 2);

    path_to_source = fileparts(mfilename('fullpath'));
    
    vol_bottom = max(0, load_mrc(fullfile(path_to_source, '../../data/FakeKvMap/FakeKvMapBeta.mrc')));
    vol_middle = max(0, load_mrc(fullfile(path_to_source, '../../data/FakeKvMap/FakeKvMapAlphaOne.mrc')));

    vol = vol_bottom + vol_middle;
    
    shifted_vols = zeros(L,L,L,n);
    for i=1:n
        shifted_fullres = fakekv_bottom_shift(vol, xyshifts(1,i), xyshifts(2,i));
        shifted_lowres = vol_downsample(shifted_fullres, L);
        shifted_vols(:,:,:,i) = shifted_lowres;
    end
end

