% SKEWED_RAND_ROTS Draw random elements from a skewed distribution over SO(3)
%
% Usage
%    rots = skewed_rand_rots(n, delta, R);
%
% Input
%    n: The number of random matrices to generate.
%    delta: The factor determining the amount of skew, where 1 gives a uniform
%       distribution, values smaller than 1 give a skew away from R, and
%       values larger than 1 give a skew towards R (default 4).
%    R: The rotation away from/toward which the distribution skews (default
%       `eye(3)`).
%
% Output
%    rots: A 3-by-3-by-n array containing as submatrices the rotation matrices
%       of n rotations sampled from the skewed distribution over SO(3) defined
%       by
%
%          R*Q(alpha, beta, gamma)
%
%      where
%
%          alpha ~ U[0, 2*pi]
%          beta ~ acos(2*U[0, 1]^delta - 1)
%          gamma = U[0, 2*pi]
%
%      are the skewed Euler angles.
%
% Note
%    The 'skewed_rand_rots' function depends on the random number state of
%    'rand', so to obtain reproducible results, its state must be controlled
%    prior to calling.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function rots = skewed_rand_rots(n, delta, R)
    if nargin < 2 || isempty(delta)
        delta = 4;
    end

    if nargin < 3 || isempty(R)
        R = eye(3);
    end

    angles = [rand(1, n)*2*pi; ...
              acos(2*rand(1, n).^delta-1); ...
              rand(1, n)*2*pi];

    rots = angles_to_rots(angles);

    rots = matfun(@(Q)(R*Q), rots, 3);
end
