function dist = geodesicSpherical(theta1,phi1,theta2,phi2, bubbleRadius)
    latitude1 = theta1-pi/2;
    longitude1 = phi1 -pi;
    latitude2 = theta2-pi/2;
    longitude2 = phi2 -pi;
    dist = geodesicll(latitude1,longitude1,latitude2,longitude2, bubbleRadius);
end
