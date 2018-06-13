data = load('startConfig3.mat');
particles = (data.centers./res)*(2*pi/360);
radii = data.radii;


sphereRadius = data.bubbleRadius;
scale = lScale(sphereRadius);
%particles = rand(10000,2);
%density = 80;
%div = 1/density;
%numStart = (density-1)*density;
%particles = zeros(numStart,2);
past = [];
%radii = ones(numStart,1)*30;
time = 130;
stepSize = 1;

radius = zeros(time/stepSize,1);

index=1;
%{
for i=div:div:1-div
    for j=0:div:1-div
        particles(index,1) = i;
        particles(index,2) = j;
        index = index+1;
    end
end

particles(:,1) = particles(:,1)*pi;
particles(:,2) = particles(:,2)*2*pi;
%}
%tic
for j=1:stepSize:time
    j
    past = [past;particles];
    %[particles,radii] =
    %step(particles,radii,scale,sphereRadius,stepSize)
    numParticles = length(particles);
    particles = gpuArray(particles);
    param1 = gpuArray(particles(:,1));
    param2 = gpuArray(particles(:,2));
    param3 = gpuArray(radii);
    param4 = gpuArray(scale*ones(numParticles,1));
    param5 = gpuArray(sphereRadius*ones(numParticles,1));
    param6 = gpuArray(stepSize*ones(numParticles,1));
    param7 = gpuArray(preloadls);
    tic
    [particlex,particley,radii] = arrayfun(@stepGPU,param1,param2,param3,param4,param5,param6,param7);
    toc
    changex = param1 - particlex;
    changey = param2 - particley;
    particles = [gather(particlex),gather(particley)];
    tic
    [particles,radii] = merge(particles,gather(radii),sphereRadius);
    toc
    %plotSphere(particles,sphereRadius);

    radius(round((j-1)/stepSize+1)) = mean(radii);
end
%{
hold off
numPast = length(past);
plotSphere(particles,sphereRadius);
hold on
%{
for i=1:numParticles
    plotIsland(past(i,1),past(i,2),radii(1),sphereRadius,scale)
end
%}
    x = sin(past(:,1)).*cos(past(:,2)).*sphereRadius;
    y = sin(past(:,1)).*sin(past(:,2)).*sphereRadius;
    z = cos(past(:,1)).*sphereRadius;
scatter3(x,y,z)
hold off
%}
res = 10;
imgSizeTheta = 180*res;
imgSizePhi = 360*res;
fImage = zeros(imgSizeTheta,imgSizePhi,'logical');
fImage=gpuArray(fImage);
[Y,X] = meshgrid(1:imgSizePhi,1:imgSizeTheta);
X = gpuArray(X);
Y = gpuArray(Y);
centers = [particles(:,1)*180*res/pi,particles(:,2)*180*res/pi];
%fImage = drawFlatBubble(fImage,X,Y,centers,radii,res,bubbleRadius);
%imshow(fImage);
%toc
plot(1:130,radius)
function [particles,radii] = step(particles,radii,scale,sphereRadius,stepSize)
    diffusion = getDiff(radii);
    numParticles = length(particles);
    direction = rand(numParticles,1)*2*pi;
    %push = cos(direction).*sqrt(4*stepSize.*diffusion)./scale;
    particles(:,1) = particles(:,1) + sin(direction).*sqrt(4*stepSize.*diffusion)./scale;
    particles(:,2) = particles(:,2) + cos(direction).*sqrt(4*stepSize.*diffusion)./scale;
    
    particles =standardize(particles);
    [particles,radii] = merge(particles,radii,sphereRadius);
    
    
    
    
end

function ls = preloadls
    layers = 20;
    length = 3.17E-3;
    hFilm = layers*length;
    visMat = 0.052*2;
    visAir = 1.827E-5;
    
    ls = hFilm*visMat/(2*visAir);
end
function [particlex,particley,radius] = stepGPU(particlex,particley,radius,scale,sphereRadius,stepSize,ls)
    
    %getDiff(radius);
    %%%%%%%%%%%%%%%
    layers = 20;
    length = 3.17E-3;
    hFilm = layers*length;
    visMat = 0.052*2;
    
    gamma = 0.577215;
    T = 303;
    kb = 1.380648E-23;
    b1 = 2.74819;
    b2 = 0.51465;
    c1 = 0.73761;
    c2 = 0.52119;
    
    
    %mob = (1./(4.*pi.*visMat.*hFilm)).*(log(ls./radii-gamma));
    %mob = 1./(16.*visAir.*radii);
    ep = radius./ls;
    mob = (1/(4*pi*visMat*hFilm))*((log(2./ep)-gamma+4.*ep./pi-(ep^2/2)*log(2/ep))/(1-(ep^3/pi)*log(2/ep)+(c1*ep^b1)/(1+c2*ep^b2)));
    %redmob = ((log(2/ep)-gamma+4*ep/pi-(ep^2/2)*log(2/ep))/(1-(ep^3/pi)*log(2/ep)+(c1*ep^b1)/(1+c2*ep^b2)));
    diffusion = mob*kb*T*(1E6)^3;
    
    
    %%%%%%%%%%%%%%%%%%%%%
    direction = rand('double')*2*pi;
    %push = cos(direction).*sqrt(4*stepSize.*diffusion)./scale;
    [particlex,particley] = reverseHaversineSpherical(particlex,particley,direction,sqrt(4*stepSize.*diffusion),sphereRadius);

    

        if(particlex > 2*pi)
            particlex = particlex-2*pi;
        elseif (particlex < 0)
            particlex = particlex+2*pi;
        end
        if(particley > 2*pi)
            particley = particley-2*pi;
        elseif (particley < 0)
            particley = particley+2*pi;
        end
    

    
    
    
    
end

function particles = standardize(particles)
    numParticles = length(particles);
    
    for i=1:numParticles
        if(particles(i,1) > 2*pi)
            particles(i,1) = particles(i,1)-2*pi;
        elseif (particles(i,1) < 0)
            particles(i,1) = particles(i,1)+2*pi;
        end
        if(particles(i,2) > 2*pi)
            particles(i,2) = particles(i,2)-2*pi;
        elseif (particles(i,2) < 0)
            particles(i,2) = particles(i,2)+2*pi;
        end
       
    end
end

function particle = standardizeGPU(particle)


        if(particle(1) > 2*pi)
            particle(1) = particle(1)-2*pi;
        elseif (particle(1) < 0)
            particle(1) = particle(1)+2*pi;
        end
        if(particle(2) > 2*pi)
            particle(2) = particle(2)-2*pi;
        elseif (particle(1) < 0)
            particle(2) = particle(2)+2*pi;
        end
       
end

function diff = getDiff(radii)
    layers = 5;
    length = 3.17E-3;
    hFilm = layers*length;
    visMat = 0.052;
    visAir = 1.827E-5;
    gamma = 0.577215;
    T = 303;
    kb = 1.380648E-23;
    b1 = 2.74819;
    b2 = 0.51465;
    c1 = 0.73761;
    c2 = 0.52119;
    
    ls = hFilm*visMat/(2*visAir);
    %mob = (1./(4.*pi.*visMat.*hFilm)).*(log(ls./radii-gamma));
    %mob = 1./(16.*visAir.*radii);
    ep = radii./ls;
    mob = (1./(4.*pi.*visMat.*hFilm)).*((log(2./ep)-gamma+4.*ep./pi-(ep.^2./2).*log(2./ep))./(1-(ep.^3./pi).*log(2./ep)+(c1.*ep.^b1)./(1+c2.*ep.^b2)));
    redmob = ((log(2./ep)-gamma+4.*ep./pi-(ep.^2./2).*log(2./ep))./(1-(ep.^3./pi).*log(2./ep)+(c1.*ep.^b1)./(1+c2.*ep.^b2)));
    diff = mob*kb*T*(1E6).^3;
end


function [particles,radii] = merge(particles,radii,sphereRadius)

    numParticles = length(particles);
    
    i = 1;
    while i<=numParticles       
        j=i+1;
        while j<=numParticles
            dist = geodesicSpherical(particles(i,1),particles(i,2),particles(j,1),particles(j,2),sphereRadius);
            minDist = radii(i)+radii(j);
            if(minDist >= dist)
                if(particles(i,2)>3*pi/2 && particles(j,2)<pi/2)
                    particles(i,1) = (particles(i,1)*radii(i)^2+particles(j,1)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                    particles(i,2) = (particles(i,2)*radii(i)^2+(particles(j,2)+2*pi)*radii(j)^2)/(radii(i)^2+radii(j)^2); 
                    if (particles(i,2) > 2*pi)
                        particles(i,2) = particles(i,2)-2*pi;
                    end
                elseif(particles(j,2)>3*pi/2 && particles(i,2)<pi/2)
                    particles(i,1) = (particles(i,1)*radii(i)^2+particles(j,1)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                    particles(i,2) = ((particles(i,2)+2*pi)*radii(i)^2+particles(j,2)*radii(j)^2)/(radii(i)^2+radii(j)^2); 
                    if (particles(i,2) > 2*pi)
                        particles(i,2) = particles(i,2)-2*pi;
                    end  
                else
                    particles(i,1) = (particles(i,1)*radii(i)^2+particles(j,1)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                    particles(i,2) = (particles(i,2)*radii(i)^2+particles(j,2)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                end
                radii(i) = sqrt(radii(i)^2+radii(j)^2);
                particles(j,:) = [];
                radii(j) = [];
                numParticles = numParticles-1;
                j=j-1;
                
            end
            j=j+1;
        end
        i=i+1;
    end

end


function [particles,radii] = mergeGPU(particles,radii,sphereRadius)

    numParticles = length(particles);
    radii = gpuArray(radii);
    particles = gpuArray(particles);
    
    
    i = 1;
    while i<=numParticles       
        j=i+1;
        dist = arrayfun(@geodesicSpherical,particles(i:end,1),particles(i:end,2),particles(i,1),particles(i,2),sphereRadius);
        minDist = radii+radii(i);
        merge = dist<=minDist;
        while j<=numParticles
            dist = geodesicSpherical(particles(i,1),particles(i,2),particles(j,1),particles(j,2),sphereRadius);
            minDist = radii(i)+radii(j);
            if(minDist >= dist)
                if(particles(i,2)>3*pi/2 && particles(j,2)<pi/2)
                    particles(i,1) = (particles(i,1)*radii(i)^2+particles(j,1)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                    particles(i,2) = (particles(i,2)*radii(i)^2+(particles(j,2)+2*pi)*radii(j)^2)/(radii(i)^2+radii(j)^2); 
                    if (particles(i,2) > 2*pi)
                        particles(i,2) = particles(i,2)-2*pi;
                    end
                elseif(particles(j,2)>3*pi/2 && particles(i,2)<pi/2)
                    particles(i,1) = (particles(i,1)*radii(i)^2+particles(j,1)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                    particles(i,2) = ((particles(i,2)+2*pi)*radii(i)^2+particles(j,2)*radii(j)^2)/(radii(i)^2+radii(j)^2); 
                    if (particles(i,2) > 2*pi)
                        particles(i,2) = particles(i,2)-2*pi;
                    end  
                else
                    particles(i,1) = (particles(i,1)*radii(i)^2+particles(j,1)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                    particles(i,2) = (particles(i,2)*radii(i)^2+particles(j,2)*radii(j)^2)/(radii(i)^2+radii(j)^2);
                end
                radii(i) = sqrt(radii(i)^2+radii(j)^2);
                particles(j,:) = [];
                radii(j) = [];
                numParticles = numParticles-1;
                j=j-1;
                
            end
            j=j+1;
        end
        i=i+1;
    end

end


function [lat,lon] = sToLL(particles)
    lat = (360/(2*pi))*(particles(:,1)-pi/2);
    lon = (360/(2*pi))*(particles(:,2) -pi);
end

function plotSphere(particles,sphereRadius)
    x = sin(particles(:,1)).*cos(particles(:,2)).*sphereRadius;
    y = sin(particles(:,1)).*sin(particles(:,2)).*sphereRadius;
    z = cos(particles(:,1)).*sphereRadius;
    
    [sx,sy,sz] = sphere;
    sx = sx*sphereRadius;
    sy = sy*sphereRadius;
    sz = sz*sphereRadius;
    
    surf(sx,sy,sz)
    hold on
    scatter3(x,y,z,'r*')
    hold off
end

function plotIsland(theta,phi,radiusIsl,sphereRadius,scale)
%{    
    xp = 0:0.001:1;
    yp = 0:0.001:1;
    [fillerp,fillert] = meshgrid(xp,yp);
    fillert = (fillert-0.5)*scale*radius-theta;
    fillerp = (fillerp-0.5)*scale*radius-phi;
    
    x = sin(fillert).*cos(fillerp).*sphereRadius;
    y = sin(fillert).*sin(fillerp).*sphereRadius;
    z = cos(fillert).*sphereRadius;
    
    scatter3(x,y,z,'r*')
    %}
    
    

    sizeC = 20;
    [columnsInImage,rowsInImage] = meshgrid(1:sizeC, 1:sizeC);
    centerX = sizeC/2;
    centerY = sizeC/2;
    radius = sizeC/2;
    circlePixels = (rowsInImage - centerY).^2 ...
        + (columnsInImage - centerX).^2 <= radius.^2;
    % circlePixels is a 2D "logical" array.
    % Now, display it.

    x=[];
    y=[];
    for i=1:sizeC
        for j=1:sizeC
            if circlePixels(i,j)
                x = [x;i];
                y = [y;j];
            end
        end
    end
    
    xgh =((x-radius)/radius)*radiusIsl/(scale*2);
    ygh = ((y-radius)/radius)*radiusIsl/scale;
    
    x = sin(theta+xgh).*cos(phi+ygh).*sphereRadius;
    y = sin(theta+xgh).*sin(phi+ygh).*sphereRadius;
    z = cos(theta+xgh).*sphereRadius;

    scatter3(x,y,z,0.1,'r*')


end


function [a,c,dlat,dlon]=haversine(lat1,lon1,lat2,lon2)
% HAVERSINE_FORMULA.AWK - converted from AWK 
    dlat = lat2-lat1;
    dlon = lon2-lon1;
    lat1 = lat1;
    lat2 = lat2;
    a = (sin(dlat./2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
    c = 2 .* asin(sqrt(a));
    %arrayfun(@(x) printf("distance: %.4f km\n",6372.8 * x), c);
end

function dist = geodesicxy(x1,y1,x2,y2, bubbleRadius, bubbleCenterx, bubbleCentery)


    x1a=x1-bubbleCentery;
    y1a=y1-bubbleCenterx;
    
    z1a = (bubbleRadius^2-x1a^2-y1a^2)^(1/2);
    longitude1 = atan2d(y1a,z1a);
    latitude1 = -asind(x1a/bubbleRadius);
    
    x2a=x2-bubbleCentery;
    y2a=y2-bubbleCenterx;
    
    z2a = (bubbleRadius^2-x2a^2-y2a^2)^(1/2);
    longitude2 = atan2d(y2a,z2a);
    latitude2 = -asind(x2a/bubbleRadius);
    
    [a,c,dlat,dlon]=haversine(latitude1,longitude1,latitude2,longitude2);
    dist = bubbleRadius*c;
end

function dist = geodesicll(latitude1,longitude1,latitude2,longitude2, bubbleRadius)

    [a,c,dlat,dlon]=haversine(latitude1,longitude1,latitude2,longitude2);
    dist = bubbleRadius*c;
end


function dist = geodesicSpherical(theta1,phi1,theta2,phi2, bubbleRadius)
    latitude1 = theta1-pi/2;
    longitude1 = phi1 -pi;
    latitude2 = theta2-pi/2;
    longitude2 = phi2 -pi;
    dist = geodesicll(latitude1,longitude1,latitude2,longitude2, bubbleRadius);
end


function scale = lScale(radius)

    circumference = 2*pi*radius;
    scale = circumference/(2*pi);

end

function rad = radians(degree) 
% degrees to radians
    rad = degree .* pi / 180;
end



