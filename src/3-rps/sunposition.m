function [x ,y, z]=sunposition(elevation,azimuth)
    x=cosd(elevation).*sind(azimuth);
    y=cosd(elevation).*cosd(azimuth);
    z=sind(elevation);

    