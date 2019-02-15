function DisplayModelPdf( fgmodel , bgmodel )

    % display foreground-background pdf
    if ~isempty( fgmodel ) && ~isempty( bgmodel )

        minMean = min( [ fgmodel.mean ; bgmodel.mean ] );
        maxMean = max( [ fgmodel.mean ; bgmodel.mean ] );
        minSigma = sqrt( min( [ fgmodel.covariance(:) ; bgmodel.covariance(:) ] ) );
        maxSigma = sqrt( max( [ fgmodel.covariance(:) ; bgmodel.covariance(:) ] ) );

        x = ( ( minMean - 3 * maxSigma ) : minSigma/20 : ( maxMean + 3 * maxSigma ) )';

        figure;

        hold on;

            plot( x , pdf( fgmodel.gmobj , x ) , 'gx-' , 'LineWidth' , 2.0 );

            plot( x , pdf( bgmodel.gmobj , x ) , 'rx-' , 'LineWidth' , 2.0 );

        hold off;

        title( 'Background and Foreground PDFs modeling using Mixture of Gaussians' , ...
               'FontSize' , 12 , 'FontWeight' , 'Bold' );

        set( gcf , 'Name' , 'Background and Foreground PDFs modeling using Mixture of Gaussians' );

        legend( { 'Foreground pdf' , 'Background pdf' } );    

    end    

end