% UNIF_RAND_ROTS Uniformly enerate random elements of SO(3)
%
% Usage
%    rots = unif_rand_rots(n);
%
% Input
%    n: The number of random matrices to generate.
%
% Output
%    rots: A 3-by-3-by-n array containing as submatrices the rotation matrices
%       of n rotations sampled uniformly over SO(3).
%
% Note
%    The 'unif_rand_rots' function depends on the random number state of
%    'rand', so to obtain reproducible results, its state must be controlled
%    prior to calling.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function rots = unif_rand_rots(n)
    angles = [rand(1, n)*2*pi; ...
              acos(2*rand(1, n)-1); ...
              rand(1, n)*2*pi];

    rots = angles_to_rots(angles);
end
