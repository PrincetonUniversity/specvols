% FB_BASIS_EVALUATE Evaluate coefficient vector in Fourier-Bessel basis
%
% Usage
%    x = fb_basis_evaluate(basis, v);
%
% Input/Output
%    See documentation for `basis_evaluate`.
%
% See also
%    basis_evaluate

% Author
%    Joakim Anden <janden@flatironinstitute.org>

function x = fb_basis_evaluate(basis, v)
    d = numel(basis.sz);

    if d == 2
        x = fb_basis_evaluate_2d(basis, v);
    elseif d == 3
        x = fb_basis_evaluate_3d(basis, v);
    end
end

function x = fb_basis_evaluate_2d(basis, v)
    [v, sz_roll] = unroll_dim(v, 2);

    r_unique = basis.coords.r_unique;
    ang_unique = basis.coords.ang_unique;
    r_idx = basis.coords.r_idx;
    ang_idx = basis.coords.ang_idx;
    mask = basis.coords.mask;

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    x = zeros([prod(basis.sz) size(v, 2)], class(v));
    for ell = 0:basis.ell_max
        k_max = basis.k_max(ell+1);

        idx_radial = ind_radial + [0:k_max-1];

        nrms = basis.norms(idx_radial);

        if ~is_precomp
            radial = zeros(size(r_unique, 1), k_max);
            for k = 1:k_max
                radial(:,k) = besselj(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        if ell == 0
            sgns = 1;
        else
            sgns = [1 -1];
        end

        for sgn = sgns
            if ~is_precomp
                if sgn == 1
                    ang = cos(ell*ang_unique);
                else
                    ang = sin(ell*ang_unique);
                end
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:k_max-1];

            x(mask,:) = x(mask,:) + ang_radial*v(idx,:);

            ind = ind + numel(idx);

            ind_ang = ind_ang + 1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end

function x = fb_basis_evaluate_3d(basis, v)
    [v, sz_roll] = unroll_dim(v, 2);

    r_unique = basis.coords.r_unique;
    ang_unique = basis.coords.ang_unique;
    r_idx = basis.coords.r_idx;
    ang_idx = basis.coords.ang_idx;
    mask = basis.coords.mask;

    is_precomp = isfield(basis, 'precomp');

    ind = 1;

    ind_radial = 1;
    ind_ang = 1;

    x = zeros([prod(basis.sz) size(v, 2)], class(v));
    for ell = 0:basis.ell_max
        idx_radial = ind_radial + [0:basis.k_max(ell+1)-1];

        nrms = basis.norms(idx_radial);

        if ~is_precomp
            radial = zeros(size(r_unique, 1), numel(idx_radial));
            for k = 1:numel(idx_radial)
                radial(:,k) = sph_bessel(ell, basis.r0(k,ell+1)*r_unique);
            end
        else
            radial = basis.precomp.radial(:,idx_radial);
        end

        radial = bsxfun(@times, radial, 1./nrms');

        for m = -ell:ell
            if ~is_precomp
                ang = real_sph_harmonic(ell, m, ang_unique(:,1), ...
                    ang_unique(:,2));
            else
                ang = basis.precomp.ang(:,ind_ang);
            end

            ang_radial = bsxfun(@times, ang(ang_idx), radial(r_idx,:));

            idx = ind + [0:numel(idx_radial)-1];

            x(mask,:) = x(mask,:) + ang_radial*v(idx,:);

            ind = ind + numel(idx);
            ind_ang = ind_ang+1;
        end

        ind_radial = ind_radial + numel(idx_radial);
    end

    x = reshape(x, [basis.sz size(x, 2)]);

    x = roll_dim(x, sz_roll);
end
