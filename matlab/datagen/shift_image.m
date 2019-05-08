function shifted_image = shift_image(image, x_shift, y_shift, extrapval)
    [nrow, ncol] = size(image);
    [X,Y] = meshgrid(1:ncol, 1:nrow);
    shifted_image = interp2(image, X-x_shift, Y-y_shift, 'linear', extrapval);
end

