function [ imSegMask_kbmaxflow ] = segObject2D_InteractiveGraphCuts( imInput , DisplayRange )
%% segObject2D_InteractiveGraphCuts allows the user to interactively segment an image into foreground and background regions
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
%     A 2D matrix containing the grayscale image that needs to be segmented.
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
%     imMask = segObject2D_maxflow();    
%     
% -- Example-1
% 
%     I = imread('coins.png');
%     seg_mask = segObject2D_maxflow( I , [0 255] );
% 
    
    if ~exist( 'imInput' , 'var' )
        
       clc
       clear all
       close all
       fclose all
       
       imInput = zeros( 100 , 100 );
       
       gmobj_bgnd_ex = gmdistribution( 50 , 30 );
       imInput_bgnd_ex = reshape( random( gmobj_bgnd_ex , numel( imInput ) ) , size( imInput ) );

       gmobj_fgnd_ex = gmdistribution( 90 , 30 );
       imInput_fgnd_ex = reshape( random( gmobj_fgnd_ex , numel( imInput ) ) , size( imInput ) );
       
       imInput = imInput_bgnd_ex;
       
       image_size = size( imInput );
       
       imInput( image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75 ) = imInput_fgnd_ex( image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75 );
       
       imInput = round( imInput );
       
    end
    
    if ~exist( 'DisplayRange' , 'var' )        
        DisplayRange = [];        
    end

    imInput = double( imInput );

    SEED_NEIGH = 3;    
    
    imsize = size( imInput );
    
    % get seed points from user
    [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_strokes( imInput, DisplayRange );

    % Build foreground and background pdf
    fgmodel = BuildModelPdf( fgnd_seed_points , imInput, SEED_NEIGH );
    bgmodel = BuildModelPdf( bgnd_seed_points , imInput , SEED_NEIGH );    

        % display foreground-background pdf
        DisplayModelPdf( fgmodel , bgmodel );
        
        % Display seed-image overlay        
        figure, imageMaskOverlay( imInput, cat( 3, fgmodel.im_seed_mask,  bgmodel.im_seed_mask ), [ 0 1 0 ; 0 0 1 ] , [ 0.7 , 0.7 ], [] );
        title( 'Seed Overlay', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Seed Overlay' );                                 
    
    % Setup N-Links
    img_ind = ( 1:prod( imsize ) )';
    [ img_ind_y , img_ind_x ] = ind2sub( imsize , img_ind );

    ndir_enode1 = [];
    ndir_enode2 = [];
    ndir_ecap = [];

        % prepare DOG kernels
        sigma = 0.5;
        gaussKernel = fspecial('gauss', round( 2 * (3 * sigma) ) * [1 1] + 1, sigma);    

        DoGKernel_x = conv2( gaussKernel , fspecial( 'sobel' )' , 'valid' );       
        DoGKernel_y = conv2( gaussKernel , fspecial( 'sobel' ) , 'valid' );

        % compute gradient-based costs
        GradX = imfilter( imInput , DoGKernel_x , 'same' , 'conv' , 'replicate' );
        GradY = imfilter( imInput , DoGKernel_y , 'same' , 'conv' , 'replicate' );

        GradXCost = exp( -GradX.^2 / ( 2 * cov( GradX( fgmodel.im_seed_mask > 0 | bgmodel.im_seed_mask > 0 ) ) ) );
        GradYCost = exp( -GradY.^2 / ( 2 * cov( GradY( fgmodel.im_seed_mask > 0 | bgmodel.im_seed_mask > 0 ) ) ) );
        
            % xdir        
            img_ind_xdir = img_ind( img_ind_x + 1 <= imsize(2) );
            img_ind_xdir_neigh = sub2ind( imsize , img_ind_y( img_ind_xdir ) , img_ind_x( img_ind_xdir ) + 1 );

            ndir_enode1 = [ ndir_enode1 ; img_ind_xdir ];            
            ndir_enode2 = [ ndir_enode2 ; img_ind_xdir_neigh ];            
            ndir_ecap   = [ ndir_ecap   ; GradXCost(img_ind_xdir) ];

            clear img_ind_xdir;
            clear img_ind_xdir_neigh;

            %ydir       
            img_ind_ydir = img_ind( img_ind_y + 1 <= imsize(1) );
            img_ind_ydir_neigh = sub2ind( imsize , img_ind_y( img_ind_ydir ) + 1 , img_ind_x( img_ind_ydir ) );

            ndir_enode1 = [ ndir_enode1 ; img_ind_ydir ];            
            ndir_enode2 = [ ndir_enode2 ; img_ind_ydir_neigh ];            
            ndir_ecap   = [ ndir_ecap   ; GradYCost(img_ind_ydir)  ];

            clear img_ind_ydir;
            clear img_ind_ydir_neigh;

            % Display Grad costs
            figure, imshow( GradXCost , [] );
            title( 'Second-Order Clique Potentials -- XDir', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
            set( gcf , 'Name' , 'Second-Order Clique Potentials -- XDir' );         

            figure, imshow( GradYCost , [] );
            title( 'Second-Order Clique Potentials -- YDir', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
            set( gcf , 'Name' , 'Second-Order Clique Potentials -- YDir' );         

        % compute second-order cliques
        SecondOrderCliques = [ ndir_enode1 , ndir_enode2 ];

        SecondOrderCliquePotential = zeros( 2 , 2 , numel( ndir_enode1 ) );

        %ndir_ecap = NormalizeMatrix( ndir_ecap );

        SecondOrderCliquePotential( 1 , 2 , : ) = ndir_ecap;
        SecondOrderCliquePotential( 2 , 1 , : ) = ndir_ecap;            

    % compute appearance term
    
        % Background
        appearance_term_bgnd = reshape( -log( pdf( bgmodel.gmobj , imInput(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxBgndVal = max( appearance_term_bgnd( ~isinf( appearance_term_bgnd ) ) );
            appearance_term_bgnd( isinf( appearance_term_bgnd ) ) = maxBgndVal + 1;

            figure, imshow( appearance_term_bgnd , [] );
            title( 'Appearance Term -- Background', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
            set( gcf , 'Name' , 'Appearance Term -- Background' );         


        % Foreground
        appearance_term_fgnd = reshape( -log( pdf( fgmodel.gmobj , imInput(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxFgndVal = max( appearance_term_fgnd( ~isinf( appearance_term_fgnd ) ) );
            appearance_term_fgnd( isinf( appearance_term_fgnd ) ) = maxFgndVal + 1;

            figure, imshow( appearance_term_fgnd , [] );
            title( 'Appearance Term -- Foreground', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
            set( gcf , 'Name' , 'Appearance Term -- Foreground' );         

    % Setup First-order clique potentials

        % Background
        FirstOrderCliquePotential_Bgnd = appearance_term_bgnd;

        % Foreground
        FirstOrderCliquePotential_Fgnd = appearance_term_fgnd;

    FirstOrderCliquePotential = [ FirstOrderCliquePotential_Bgnd(:) , FirstOrderCliquePotential_Fgnd(:) ];
            
    % Solve Maxflow using K-B maxflow wrapper
    fprintf( 1 , '\n\nSolving the Maximum Flow Problem using K-B Maxflow wrapper ...\n\n'  );

        tic

        [ energy , seg_labels ] = MinimizeBinaryOrder2MRFEnergy_GraphCuts( FirstOrderCliquePotential, SecondOrderCliques, SecondOrderCliquePotential );       

        energy

        imSegMask_kbmaxflow = reshape( seg_labels , imsize );
        imSegMask_kbmaxflow = imSegMask_kbmaxflow > 0;

        toc   
   
       % Display Segmentation Results
       figure; 

       imageMaskOverlay( imInput , cat( 3 , imSegMask_kbmaxflow, fgmodel.im_seed_mask > 0 , bgmodel.im_seed_mask > 0 ) , [ 0 1 0 ; 0 0 1 ; 1 0 0 ] , [ 0.3 , 0.5 , 0.5 ] , DisplayRange );

       title( 'Final Segmentation Result with seed regions overlayed (K-B Maxflow)' , ...
              'FontSize' , 12 , 'FontWeight' , 'Bold' );                       
        
end