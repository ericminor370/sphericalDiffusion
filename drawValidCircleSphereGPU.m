function fImage = drawValidCircleSphereGPU(fImagePix,x,y,seedx,seedy,radius,res,bubbleRadius)
    fImage = fImagePix&&(geodesicSphericalDegRes(seedx,seedy,x,y,bubbleRadius,res)>=radius);
    %fImage = fImagePix&&(sqrt((x-seedx)^2+(y-seedy)^2)>=radius);
end