function [ imSegMask_kbmaxflow ] = segObject3D_InteractiveGraphCuts( imInput , DisplayRange )
%% segObject3D_InteractiveGraphCuts allows the user to interactively segment an image into foreground and background regions
% This function implements the traditional interactive image segmentation using graph-cuts (See Ref[1]). 
% 
% As soon as the function is called, a GUI window will popup asking the user to interactively provide the foreground and background seed regions via 2D strokes. 
% The strokes provided using this GUI are used to build appearance models for the foreground and the background object. 
% Both foreground and background appearance is modeled using a mixture of gaussians. Also, one gaussian/normal pdf is constructed for each stroke provided by the user.
% 
% After obtaining the seed regions, a 2-label MRF problem is setup and solved to obtain the segmentation.
% 
% Sites of the MRF are pixels in the image
% Labels - {foregound, background}
% Cliques - We only consider first- and second-order cliques. The second-order cliques are defined based on a 4-neighborhood graph. 
%           So, there are two types of order-2 cliques: (i) between neighboring pixels in x-direction, and (ii) between neighboring pixels in y-direction
% 
% The first-order clique potentials corresponding to foreground and background are determined using the mixture of guassian models built from user-define seed regions. 
% 
% V_i( L_i = 1 ) = -log( pdf( I(i) , mu_fgnd ; var_fgnd ) );
% V_i( L_i = 0 ) = -log( pdf( I(i) , mu_bgnd ; var_bgnd ) );
% 
% The second-order clique potentials are defined as a generalized potts interaction model based on a function that is inversely proportional to the 
% intensity gradient between the two pixels involved in the clique.
% 
% V_ij( L_i , L_j ) = 0 if L_i = L_j
%                   = exp( -( I(i) - I(j) )^2 / (2 * sigma^2) ) if L_i ~= L_j 
% 
% As might be evident, the above defined second-order clique is submodular and hence can be minimized using graph-cuts.
% We use the implementation of Boykov and Komogrov to solve the maxflow problem (http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html).
% 
% Refer to the following papers for more details on the underlying the theory behind this implementation:
% 
% [1] Y. Boykov and G. Funka-Lea. Graph cuts and efficient ND image segmentation. International Journal of Computer Vision, 70(2):109–131, 2006.
% [2] V. Kolmogorov and R. Zabin. What energy functions can be minimized via graph cuts? IEEE Trans. on Pattern Analysis and Machine Intelligence, 26(2):147–159, 2004.
% [3] Y. Boykov, V. Kolmogorov. An experimental comparison of min-cut/max-flow algorithms for energy minimization in vision. IEEE Trans Pattern Analysis and Machine Intelligence, 26:1124–1137, 2004
% 
% Input arguments:
% 
% --  imInput: 
% 
%     A 3D matrix containing the grayscale image that needs to be segmented.
%     
%     If imInput is not provided then segmentation is performed on a synthetically generated image. This gives you a demo of what the function does and how it works.
%     
% -- DisplayRange (optional): 
% 
%    An intensity range for better visualization of the input image. 
%    Default value [min(imInput(:) max(imInput(:))]   
%    
% Output arguments:
% 
% --  imSegMask_kbmaxflow
% 
%     A binary segmentation mask
% 
% Example usages:
% 
% -- Demo Mode
% 
%     imMask = segObject3D_InteractiveGraphcuts();    
%     
    
    if ~exist( 'imInput' , 'var' )
        
       clc
       clear all
       close all
       fclose all
       
       imInput = zeros( 100 , 100, 20 );
       
       gmobj_bgnd_ex = gmdistribution( 50 , 30 );
       imInput_bgnd_ex = reshape( random( gmobj_bgnd_ex , numel( imInput ) ) , size( imInput ) );

       gmobj_fgnd_ex = gmdistribution( 90 , 30 );
       imInput_fgnd_ex = reshape( random( gmobj_fgnd_ex , numel( imInput ) ) , size( imInput ) );
       
       imInput = imInput_bgnd_ex;
       
       image_size = size( imInput );
       
       fgBox = { image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75, ':' };
       imInput(fgBox{:}) = imInput_fgnd_ex(fgBox{:});       
       
       imInput = round( imInput );
       
    end
    
    if ~exist( 'DisplayRange' , 'var' )        
        DisplayRange = [];        
    end

    SEED_NEIGH = 3;    
    regularizationWeight = 1.0;
    spacing = ones(1,3);
    
    imInput = double( imInput );   
    imsize = size( imInput );
    
    % get seed points from user
    [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_3d_strokes( imInput, DisplayRange );

    % Build foreground and background pdf
    fgmodel = BuildModelPdf( fgnd_seed_points , imInput, SEED_NEIGH );
    bgmodel = BuildModelPdf( bgnd_seed_points , imInput , SEED_NEIGH );    

        % display foreground-background pdf
        DisplayModelPdf( fgmodel , bgmodel );
        
        % Display seed-image overlay        
        imseriesmaskshow( imInput, {fgmodel.im_seed_mask, bgmodel.im_seed_mask} );
        title( 'Seed Overlay', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Seed Overlay' );                                 
    
    % Setup unary clique potentials
    
        % Background
        FirstOrderCliquePotential_Bgnd = reshape( -log( pdf( bgmodel.gmobj , imInput(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxBgndVal = max( FirstOrderCliquePotential_Bgnd( ~isinf( FirstOrderCliquePotential_Bgnd ) ) );
            FirstOrderCliquePotential_Bgnd( isinf( FirstOrderCliquePotential_Bgnd ) ) = maxBgndVal + 1;

            imseriesshow( FirstOrderCliquePotential_Bgnd );
            title( 'Unary Clique Potential -- Background', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
            set( gcf , 'Name' , 'Unary Clique Potential -- Background' );         
            
    
        % Foreground
        FirstOrderCliquePotential_Fgnd = reshape( -log( pdf( fgmodel.gmobj , imInput(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxFgndVal = max( FirstOrderCliquePotential_Fgnd( ~isinf( FirstOrderCliquePotential_Fgnd ) ) );
            FirstOrderCliquePotential_Fgnd( isinf( FirstOrderCliquePotential_Fgnd ) ) = maxFgndVal + 1;

            imseriesshow( FirstOrderCliquePotential_Fgnd );
            title( 'Unary Clique Potential -- Foreground', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
            set( gcf , 'Name' , 'Unary Clique Potential -- Foreground' );         
            
    FirstOrderCliquePotential = [ FirstOrderCliquePotential_Bgnd(:) , FirstOrderCliquePotential_Fgnd(:) ];
    
    % Setup pairwise clique potentials
    directionLabels = { 'Y', 'X', 'Z' };
    ndir_enode1 = [];
    ndir_enode2 = [];
    ndir_ecap = [];
        
        % setup pairwise cliques along each dimension and the associate clique potentials        
        img_pixind = (1:numel(imInput))';
        img_pixsubind = cell(1,ndims(imInput));
        [img_pixsubind{:}] = ind2sub( imsize, img_pixind );
        
        for i = 1:ndims(imInput)
            
            % edge cost
            imGaussGrad = filterGaussGradND(imInput, 1.0, i, 'spacing', spacing );
            sigmaGrad = std( imGaussGrad(~fgmodel.im_seed_mask) );            
            imGradCost = exp( -imGaussGrad.^2 / (2 * sigmaGrad^2) );
            
            % from node indices
            from_node_ind = img_pixind( (img_pixsubind{i} + 1) <= imsize(i) );
            
            % to node indices
            img_neighpixsubind = cell(1,ndims(imInput));
            for j = 1:ndims(imInput)
                img_neighpixsubind{j} = img_pixsubind{j}( from_node_ind );
            end
            img_neighpixsubind{i} = img_neighpixsubind{i} + 1;
            
            to_node_ind = sub2ind( imsize, img_neighpixsubind{:} );
            
            % add to list of pairwise cliques
            ndir_enode1 = [ ndir_enode1 ; from_node_ind ];
            ndir_enode2 = [ ndir_enode2 ; to_node_ind ];
            ndir_ecap   = [ ndir_ecap   ; imGradCost(from_node_ind) ];

            imseriesshow( imGradCost );
            set( gcf , 'Name' , sprintf( 'Second-Order Clique Potentials -- %s-Dir', directionLabels{i} ) );         
            
        end
        
        SecondOrderCliques = [ ndir_enode1 , ndir_enode2 ];
        SecondOrderCliquePotential = zeros( 2 , 2 , numel( ndir_enode1 ) );
        SecondOrderCliquePotential( 1 , 2 , : ) = ndir_ecap;
        SecondOrderCliquePotential( 2 , 1 , : ) = ndir_ecap;            
        SecondOrderCliquePotential = regularizationWeight * SecondOrderCliquePotential;
        
        clear img_pixind;
        clear img_pixsubind;
    
    % Solve Maxflow using K-B maxflow wrapper
    fprintf( 1 , '\n\nSolving the Maximum Flow Problem using K-B Maxflow wrapper ... '  );

        tic
        [ energy , seg_labels ] = MinimizeBinaryOrder2MRFEnergy_GraphCuts( FirstOrderCliquePotential, SecondOrderCliques, SecondOrderCliquePotential );       
        timeElapsed = toc;

        fprintf( 1 , 'It took %.2f seconds\n\n', timeElapsed  );
        
        imSegMask_kbmaxflow = reshape( seg_labels , imsize );
        imSegMask_kbmaxflow = imSegMask_kbmaxflow > 0;
    
       % Display
       imseriesmaskshow( imInput , {imSegMask_kbmaxflow, fgmodel.im_seed_mask, bgmodel.im_seed_mask} );

       title( 'Final Segmentation Result with seed regions overlayed (K-B Maxflow)' , ...
              'FontSize' , 12 , 'FontWeight' , 'Bold' );    
       set( gcf , 'Name' , 'Final Segmentation Result with seed regions overlayed (K-B Maxflow)' );
              
end