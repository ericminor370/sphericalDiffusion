function [lat2,lon2] = reverseHaversine(lat1,lon1, bearing, dist,bubbleRadius)
    lat2 = asin(sin(lat1)*cos(dist/bubbleRadius)+cos(lat1)*sin(dist/bubbleRadius)*cos(bearing));
    lon2 = lon1 + atan2(sin(bearing)*sin(dist/bubbleRadius)*cos(lat1),cos(dist/bubbleRadius)-sin(lat1)*sin(lat2));
end