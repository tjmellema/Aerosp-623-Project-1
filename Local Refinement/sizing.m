function h = sizing(X_pt, X_BE_spline_top, X_BE_spline_bottom)
        %project the point on the wall
        [d_top, xb_top, ~] = projection(X_pt, X_BE_spline_top);
        [d_bottom, xb_bottom, ~] = projection(X_pt, X_BE_spline_top);
        [d, xb]


        
        %sizing formula
        a = 1; %sixing factor
        h_max = 2; %max cell size
        % cosine spacing between leading and trailing edge
        % linear in distnace from the airfoil
        h = a.* min(h_max, (1 + cos(xb./19)) ./ 2 + d);
end