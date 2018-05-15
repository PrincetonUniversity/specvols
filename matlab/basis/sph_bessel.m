% SPH_BESSEL Compute spherical Bessel function values
%
% Usage
%    j = sph_bessel(ell, r)
%
% Input
%    ell: The order of the spherical Bessel function.
%    r: The coordinates where the function is to be evaluated.
%
% Output
%    j: The value of j_ell at r.

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function j = sph_bessel(ell, r)
    j = zeros(size(r));

    if ell == 0
        j(r==0) = 1;
    else
        j(r==0) = 0;
    end

    j(r~=0) = sqrt(pi./(2*r(r~=0))).*besselj(ell+1/2, r(r~=0));
end
