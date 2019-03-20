function volume = gen_clock_volume(L, clock_hand_angle_rad)
    assert(mod(L,2) == 1);
    H = round((L-1)/2); % Half of the side length
    
    arrow_coords = [-0.05 0; -0.05 +0.8; +0.35 +0.4; +0.05 0.55; +0.05 0; -0.05 0];

    R = [cos(clock_hand_angle_rad) -sin(clock_hand_angle_rad);
        sin(clock_hand_angle_rad) cos(clock_hand_angle_rad)];
    rotated_arrow_coords = (R*(arrow_coords'))';
    arrow_mask2d = poly2mask(rotated_arrow_coords(:,1)*H+H+1, ...
        rotated_arrow_coords(:,2)*(-H)+H+1, 2*H+1, 2*H+1);
    
    [X,Y] = meshgrid((-H:H)/H,(-H:H)/H);

    arrow3d = repmat((arrow_mask2d),1,1,round(0.2*H));

    clockface2d = X.^2 + Y.^2 <= (0.8^2+0.1^2);
    clockstand_coords = [0 0; 0.5 -1; -0.5 -1; 0 0];
    clockstand2d = poly2mask(clockstand_coords(:,1)*(-H)+H+1, ...
        clockstand_coords(:,2)*(-H)+H+1, 2*H+1, 2*H+1); 
    clock2d = clockface2d|clockstand2d;

    hour_indicators = zeros(size(X));
    for k=0:11
        angle = 2*pi*k/12;
        hour_indicators = hour_indicators + ...
            ((X-0.65*cos(angle)).^2 + (Y-0.65*sin(angle)).^2 < 0.06^2);
    end

    clock3d = repmat(clock2d,1,1,round(0.4*H));
    hour_indicators3d = repmat(hour_indicators,1,1,round(0.05*H));
    volume = cat(3, arrow3d, hour_indicators3d, clock3d);
end
