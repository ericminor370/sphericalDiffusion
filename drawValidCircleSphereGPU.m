function fImage = drawValidCircleSphereGPU(fImagePix,x,y,seedx,seedy,radius,res,bubbleRadius)
    fImage = fImagePix&&(geodesicSphericalDegRes(seedx,seedy,x,y,bubbleRadius,res)>=radius);
end