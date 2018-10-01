% PROJBACKPROJ_VOLUME Fast FFT-based tomographic projection of a 3D volume
%
% Input
%       vol: NxNxN volume
%       direction: Projection direction vector [x,y,z]
%
% Output
%       projected_volume: NxNxN input volume after tomographic projection
%
% Description
%       This function pads the input and performs Fourier convolution
%       using the kernel obtained from projbackproj_kernel().
%
function projected_volume = projbackproj_volume(vol, direction)
    assert(all(size(direction) == [1 3]));

    [N_x, N_y, N_z] = size(vol);
    assert(N_x == N_y);
    assert(N_y == N_z);
    N = N_x;

    line_kernel_dft = projbackproj_kernel_dft(N, direction);
    assert(all(size(line_kernel_dft) == [2*N, 2*N, 2*N]));
    
    padded_vol = zeros([2*N 2*N 2*N]);
    padded_vol(1:N, 1:N, 1:N) = vol;
    
    projected_padded_vol_dft = fft3(padded_vol).*line_kernel_dft;
    projected_padded_vol = real(ifft3(projected_padded_vol_dft));
    projected_volume = projected_padded_vol(1:N, 1:N, 1:N);
end
