function dist = geodesicSphericalDegRes(theta1,phi1,theta2,phi2, bubbleRadius,res)
    theta1 = radians(theta1/res);
    phi1 = radians(phi1/res);
    theta2 = radians(theta2/res);
    phi2 = radians(phi2/res);
    
    dist = geodesicSpherical(theta1,phi1,theta2,phi2, bubbleRadius);
end