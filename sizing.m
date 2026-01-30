function h = sizing(X_pt, X_BE_spline)
        %project the point on the wall
        [d, xb, ~] = projection(X_pt, X_BE_spline);
        
        %sizing formula
        a = 1; %sizing factor
        h_max = 2; %max cell size
        % cosine spacing between leading and trailing edge
        % linear in distnace from the airfoil
        h = a.* min(h_max, (1 + cos(xb./19)) ./ 2 + d);
end