function fImage = drawCircleSphere(fImage,imgSizeTheta,imgSizePhi,seed,radius,res,bubbleRadius)

       for i=1:imgSizeTheta
            for j=1:imgSizePhi
                if geodesicSphericalDegRes(seed(1),seed(2),i,j,bubbleRadius,res)<radius
                    fImage(i,j) = 1;
                end
            end
        end

end