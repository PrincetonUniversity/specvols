% PROJBACKPROJ_KERNEL_DFT Direct computation of the DFT of a tomographic projection kernel 
%
% Input
%       N: side length of the volume
%       projection_direction_vector: Direction vector [x,y,z]
%
% Output
%       kernel_dft: 2Nx2Nx2N DFT of a tomographic projection kernel.
%
% Description
%       We take the kernel to be the line segment l(t) = vt where t goes
%       from -R to +R, and compute its discrete Fourier frequencies.
%
%       R is chosen to the largest possible such that
%       v*R stays within an NxNxN volume.
%       Hence, R = N/max(abs(normalized(v)))
%
%       F(k) = \int_{-R}^{+R} exp(-2 pi i <k,vt>/2N dt
%            = sin(pi <k,v> R/N)/(pi <k,v>)
%
function kernel_dft  = projbackproj_kernel_dft(N, projection_direction_vector)
    v = projection_direction_vector / norm(projection_direction_vector);

    % Build 2Nx2Nx2N FFT frequency indices
    % of the form [0, 1, ..., N-1, -N -N+1, ..., -1]
    [X,Y,Z] = meshgrid(0:(2*N-1), 0:(2*N-1), 0:(2*N-1));
    X = fftshift(X-N);
    Y = fftshift(Y-N);
    Z = fftshift(Z-N);

    k_dot_v = X*v(1) + Y*v(2) + Z*v(3);
    
    % Note that sinc(t) = sin(pi*t)/(pi*t)
    inv_max_abs_v = 1/max(abs(v));
    kernel_dft = sinc(k_dot_v*inv_max_abs_v)*inv_max_abs_v;
end
