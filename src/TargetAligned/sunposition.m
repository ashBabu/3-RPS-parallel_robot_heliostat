function [x, y ,z]=sunposition(alpha,A)
    x=cosd(alpha).*sind(A);
    y=cosd(alpha).*cosd(A);
    z=sind(alpha);

    