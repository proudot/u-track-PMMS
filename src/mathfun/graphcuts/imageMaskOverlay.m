function imageMaskOverlay( im , maskImages , maskColorMap , maskAlpha , ImageIntensityDisplayRange )

    numMasks  = size( maskImages , 3 );
    numColors = size( maskColorMap , 1 );
    numAlpha  = length( maskAlpha );
    
    if length( unique( [ numMasks numColors numAlpha ] ) ) ~= 1
        
       error('error: number of mask images, alpha values and colors are not the same');
        
    end
    
    if exist( 'ImageIntensityDisplayRange' , 'var' )

        imshow( im , ImageIntensityDisplayRange );
	
    else
    
        imshow( im , [] );

    end
    
    maskImages = double( maskImages );
    im = double( im );
    
    hold on;
        
    imsize = size( maskImages( : , : , 1 ) );
    
    for i = 1:numMasks
        
        imRGBCurMask = zeros( [ imsize 3 ] );
        
        for c = 1:3
            
            imRGBCurMask( : , : , c ) = maskImages( : , : , i ) * maskColorMap( i , c );
            
        end
        
        imCurAlpha = maskImages( : , : , i ) * maskAlpha( i );

        image( imRGBCurMask , 'AlphaData' , imCurAlpha );
        
        hold on;
        
    end
    
    hold off;
    
    drawnow;
    
end