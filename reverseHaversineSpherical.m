function [theta2,phi2] = reverseHaversineSpherical(theta1,phi1,bearing,dist,bubbleRadius)
    lat = theta1-pi/2;
    lon = phi1-pi;
    [lat2,lon2] = reverseHaversine(lat,lon,bearing,dist,bubbleRadius);
    theta2 = lat2+pi/2;
    phi2 = lon2+pi;  
end