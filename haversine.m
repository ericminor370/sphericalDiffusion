function [a,c,dlat,dlon]=haversine(lat1,lon1,lat2,lon2)
% HAVERSINE_FORMULA.AWK - converted from AWK 
    dlat = lat2-lat1;
    dlon = lon2-lon1;
    lat1 = lat1;
    lat2 = lat2;
    a = (sin(dlat./2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
    c = 2 .* asin(sqrt(abs(a)));
    %arrayfun(@(x) printf("distance: %.4f km\n",6372.8 * x), c);
end
