function [lat2,lon2] = reverseHaversine2(lat1,lon1,bearing,dist,bubbleRadius)
    lat2 = asin(sin(lat1)*cos(dist/bubbleRadius)+cos(lat1)*sin(dist/bubbleRadius)*cos(bearing));
    if (cos(lat2)==0)
        lon2 = lon1;
    else
        lon2 = mod(lon1-asin(sin(bearing)*sin(dist/bubbleRadius)./cos(lat2))+pi,2*pi)-pi;
    end
end