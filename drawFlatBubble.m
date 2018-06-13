function fImage = drawFlatBubble(fImage,X,Y,centers,radii,res,bubbleRadius)
    
    numIslands = length(radii);
    for i=1:numIslands
        fImage = arrayfun(@drawCircleSphereGPU,fImage,X,Y,centers(i,1),centers(i,2),radii(i),res,bubbleRadius);
    end
end