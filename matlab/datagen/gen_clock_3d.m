function volume = gen_clock_3d(L, clock_hand_angle_rad)
    assert(mod(L,2) == 1);
    H = round((L-1)/2); % Half of the side length
    
    % Render clock hand
    clock_hand_2d = gen_clock_hand_2d(L, clock_hand_angle_rad);
    clock_hand_3d = repmat((clock_hand_2d),1,1,round(0.2*H));
    
    [X,Y] = meshgrid((-H:H)/H,(-H:H)/H);

    % Render clock body
    clockface2d = X.^2 + Y.^2 <= (0.8^2+0.1^2);
    clockstand_coords = [0 0; 0.5 -1; -0.5 -1; 0 0];
    clockstand2d = poly2mask(clockstand_coords(:,1)*(-H)+H+1, ...
        clockstand_coords(:,2)*(-H)+H+1, 2*H+1, 2*H+1); 
    clock2d = clockface2d|clockstand2d;
    clock3d = repmat(clock2d,1,1,round(0.4*H));
    
    % Render hour tick marks 
    hour_indicators = zeros(size(X));
    for k=0:11
        angle = 2*pi*k/12;
        hour_indicators = hour_indicators + ...
            ((X-0.65*cos(angle)).^2 + (Y-0.65*sin(angle)).^2 < 0.06^2);
    end
    hour_indicators3d = repmat(hour_indicators,1,1,round(0.05*H));

    % Join all the components to form the complete volume
    volume = cat(3, clock3d, hour_indicators3d, clock_hand_3d);
    
    % Pad 3rd axis symmetrically
    padding_size = (L - size(volume,3));
    volume = padarray(volume, [0 0 floor(padding_size/2)], 'both');
    % Pad any leftovers non-symmetrically
    volume = padarray(volume, [0 0 L-size(volume,3)], 'pre');
end
