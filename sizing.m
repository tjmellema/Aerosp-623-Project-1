function h = sizing(X_pt, X_BE_spline_top, X_BE_spline_bottom)
        %project the point on the wall
        [d_top, xb_top, ~] = projection(X_pt, X_BE_spline_top);
        [d_bottom, xb_bottom, ~] = projection(X_pt, X_BE_spline_bottom);
        if d_top < d_bottom
            d = d_top; 
            xb = xb_top; 
        else
            d = d_bottom; 
            xb = xb_bottom; 
        end
        
        %sizing formula
        a = 1; %sizing factor
        h_max = 2; %max cell size
        h_min = 0.0001; %min cell size
        % cosine spacing between leading and trailing edge
        % linear in distnace from the airfoil
        h = a.* max(min(h_max, (1 + cos(xb./19)) ./ 2 + .25*d), h_min);
end