function dist = geodesicll(latitude1,longitude1,latitude2,longitude2, bubbleRadius)

    [a,c,dlat,dlon]=haversine(latitude1,longitude1,latitude2,longitude2);
    dist = bubbleRadius*c;
end