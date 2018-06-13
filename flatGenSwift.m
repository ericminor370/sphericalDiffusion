% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.

imgSize = 1024;
fImage = zeros(imgSize,imgSize,'logical');
radius = 100;
centers=[];
radii=[];

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
while fail~=1
    [fImage,centers,radii,fail] = addIsland(fImage,100,centers,radii,imgSize);
    %[fImage,centers,radii,fail] = addIsland(fImage,100,centers,radii,imgSize);
end
fail = 0;
while fail~=1
    [fImage,centers,radii,fail] = addIsland(fImage,50,centers,radii,imgSize);
    %[fImage,centers,radii,fail] = addIsland(fImage,100,centers,radii,imgSize);
end
fail = 0;
while fail~=1
    [fImage,centers,radii,fail] = addIsland(fImage,20,centers,radii,imgSize);
    %[fImage,centers,radii,fail] = addIsland(fImage,100,centers,radii,imgSize);
end
fail = 0;
while fail~=1
    [fImage,centers,radii,fail] = addIsland(fImage,10,centers,radii,imgSize);
    %[fImage,centers,radii,fail] = addIsland(fImage,100,centers,radii,imgSize);
end



imshow(fImage)

function [fImage,centers,radii,fail] = addIsland(fImage,radius,centers,radii,imgSize)
    valid = ones(imgSize,imgSize,'logical');
    
    for i=1:imgSize
        for j=1:imgSize
            if i<radius || j<radius || imgSize-i<radius || imgSize-j<radius
                valid(i,j) = 0;
            end
        end
    end
    
    numCircles = length(radii);
    for i=1:numCircles
        valid = ~drawCircle(~valid,imgSize,centers(i,:),radii(i)+radius+1);
    end
    %imshow(valid)
    numValid = sum(sum(valid));
    if numValid>0
        %imshow(valid)
        validPixels = zeros(numValid,2);
        index = 1;
        for i=1:imgSize
            for j=1:imgSize
                if valid(i,j) == 1
                    validPixels(index,1) = i;
                    validPixels(index,2) = j;
                    index=index+1;
                end
            end
        end

        seed = datasample(validPixels,1,1);


        fImage = drawCircle(fImage,imgSize,seed,radius);
        centers = [centers;seed];
        radii = [radii;radius];
        fail = 0;
    else
        fail = 1;
    end

end

function fImage = drawCircle(fImage,imgSize,seed,radius)

       for i=1:imgSize
            for j=1:imgSize
                if sqrt((i-seed(1)).^2+(j-seed(2)).^2)<radius
                    fImage(i,j) = 1;
                end
            end
        end

end

function in = boundarySquare(radius,imgSize)

end