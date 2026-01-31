function [Iflagged, Bflagged] = FlagEdge(node_data, IE, BE, X_BE_spline_top, X_BE_spline_bot)

Iflagged = false(size(IE,1),1);
Bflagged = false(size(BE,1),1);

% Interior edges
for i = 1:size(IE,1)

    n1 = IE(i,1);
    n2 = IE(i,2);

    Xpt1 = node_data(n1,:);
    Xpt2 = node_data(n2,:);

    h_current = norm(Xpt2 - Xpt1);
    X_midpt = 0.5*(Xpt1+Xpt2);

    h_desired = sizing(X_midpt, X_BE_spline_top, X_BE_spline_bot);

    if h_desired < h_current
        Iflagged(i) = true;
    end
end

% Boundary edges
for i = 1:size(BE,1)

    n1 = BE(i,1);
    n2 = BE(i,2);

    Xpt1 = node_data(n1,:);
    Xpt2 = node_data(n2,:);

    h_current = norm(Xpt2 - Xpt1);
    X_midpt = 0.5*(Xpt1+Xpt2);

    h_desired = sizing(X_midpt, X_BE_spline_top, X_BE_spline_bot);

    if h_desired < h_current
        Bflagged(i) = true;
    end
end

end