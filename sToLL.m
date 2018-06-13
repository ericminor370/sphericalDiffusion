function [lat,lon] = sToLL(particles)
    lat = (360/(2*pi))*(particles(:,1)-pi/2);
    lon = (360/(2*pi))*(particles(:,2) -pi);
end
