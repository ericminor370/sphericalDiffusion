% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
%poolobj = gcp;
%addAttachedFiles(poolobj,{'drawCircleSphereGPU.m','haversine.m','radians.m','geodesicll.m','geodesicSpherical.m','geodesicSphericalDegRes.m','sToLL'})
startData = load('startData');
distRadii = sort(startData.radii);
mult = 1/startData.sFract;
userMult = 1.2; %additional multiplier to account for non detected islands
%multFract = 1/startData.sFract - mult;
%numIslands = floor(lengt(distRadii*mult));
res = 10;
imgSizeTheta = 180*res;
imgSizePhi = 360*res;
fImage = zeros(imgSizeTheta,imgSizePhi,'logical');
bubbleRadius = startData.bubbleRadius;
centers=[];
radii=[];
[N,edges] = histcounts(distRadii);

binCents = (edges(1:end-1)+edges(2:end))/2;
binCents = binCents(end:-1:1);
N = floor(N(end:-1:1)*mult*userMult);
numBins = length(binCents);

%imageSizeX = radius*2;
%imageSizeY = radius*2;
%[columnsInImage rowsInImage] = meshgrid(1:imgSize, 1:imgSize);
% Next create the circle in the image.
%centerX = radius;
%centerY = radius;

%circlePixels = (rowsInImage - centerY).^2 ...
%    + (columnsInImage - centerX).^2 <= radius.^2;
% circlePixels is a 2D "logical" array.
% Now, display it.
%imshow(circlePixels) ;
%colormap([0 0 0; 1 1 1]);
%title('Binary image of a circle');
fail = 0;
%while fail~=1
tic
fImage=gpuArray(fImage);
[Y,X] = meshgrid(1:imgSizePhi,1:imgSizeTheta);
X = gpuArray(X);
Y = gpuArray(Y);
%numBins=11;
validStack = gpuArray(ones(imgSizeTheta,imgSizePhi,numBins,'logical'));
counter = 0;
failj = 0;
faili = 0;
first = 1;
totalIslands = sum(N);

for j=1:numBins
    for i=1:N(j)
        %imshow(fImage)
        [centers,radii,validStack,fail] = addIsland(j,centers,radii,imgSizeTheta,imgSizePhi,numBins,res,bubbleRadius,X,Y,validStack,binCents);
        counter = counter+1;
        counter/totalIslands
        if fail&&first
            failj = j;
            faili = i;
            first = 0;
        end
    end
end

fImage = drawFlatBubble(fImage,X,Y,centers,radii,res,bubbleRadius);

toc


imshow(fImage)
i = 1;
filename = strcat('startConfig',num2str(i),'.mat');
while exist(filename,'file') == 2
    i = i+1;
    filename = strcat('startConfig',num2str(i),'.mat');
end
save(filename,'fImage','centers','radii','userMult','bubbleRadius');



function [centers,radii,validStack,fail] = addIsland(binPos,centers,radii,imgSizeTheta,imgSizePhi,numBins,res,bubbleRadius,X,Y,validStack,binCents)
    %valid = gpuArray(ones(imgSizeTheta,imgSizePhi,'logical'));
    
    valid = gather(validStack(:,:,binPos));
    numValid = sum(sum(valid));
    if numValid>0
        %imshow(valid)
        validPixels = zeros(numValid,2);
        index = 1;
        for i=1:imgSizeTheta
            for j=1:imgSizePhi
                if valid(i,j) == 1
                    validPixels(index,1) = i;
                    validPixels(index,2) = j;
                    index=index+1;
                end
            end
        end

        seed = datasample(validPixels,1,1);
        
        for i=1:numBins
            validStack(:,:,i) = arrayfun(@drawValidCircleSphereGPU,validStack(:,:,i),X,Y,seed(1),seed(2),binCents(binPos)+binCents(i)+1,res,bubbleRadius);
        end
        
        

        %fImage = arrayfun(@drawCircleSphereGPU,fImage,X,Y,seed(1),seed(2),radius,res,bubbleRadius);
        centers = [centers;seed];
        radii = [radii;binCents(binPos)];
        fail = 0;
    else
        fail = 1;
    end

end






