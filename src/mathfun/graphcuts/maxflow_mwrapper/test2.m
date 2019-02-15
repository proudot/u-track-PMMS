
% A simple test of foreground/background segmentation using a simple binary image

clc
clear all
close all
fclose all

imInput = zeros( 512 , 512 );
imsize = size( imInput );
imInput( round( imsize(1)/4 ) : round( 3 * imsize(1) / 4 ) , round( imsize(2)/4 ) : round( 3 * imsize(2) / 4 ) ) = 255;

figure, imshow( imInput , [] );
title( 'Input Image' );
set( gcf , 'Name' , 'Input Image' );         


% Setup N-Links based on Generalized Potts Model

    % Compute directional image gradient for N-Links
    W = 1.0;

    gaussKernel = fspecial('gauss', [13 13], sqrt(13));

    DoGKernel_x = conv2( gaussKernel , fspecial( 'sobel' )' , 'valid' );       
    DoGKernel_y = conv2( gaussKernel , fspecial( 'sobel' ) , 'valid' );
    
    GradX = W * imfilter( imInput , DoGKernel_x , 'same' , 'conv' , 'replicate' );
    GradY = W * imfilter( imInput , DoGKernel_y , 'same' , 'conv' , 'replicate' );
    
    VarNLink = 0.5;
    
    % create NeighGraph
    img_ind = ( 1:prod( imsize ) )';
    [ img_ind_y , img_ind_x ] = ind2sub( imsize , img_ind );
    
    ndir_enode1 = [];
    ndir_enode2 = [];
    ndir_ecap = [];
    
        % xdir
        img_ind_xdir = img_ind( img_ind_x + 1 <= imsize(2) );
        img_ind_xdir_neigh = sub2ind( imsize , img_ind_y( img_ind_xdir ) , img_ind_x( img_ind_xdir ) + 1 );

        ndir_enode1 = [ ndir_enode1 ; img_ind_xdir ];            
        ndir_enode2 = [ ndir_enode2 ; img_ind_xdir_neigh ];            
        ndir_ecap   = [ ndir_ecap   ; exp( -GradX(img_ind_xdir).^2 / VarNLink ) ];

        clear img_ind_xdir;
        clear img_ind_xdir_neigh;
        
        %ydir
        img_ind_ydir = img_ind( img_ind_y + 1 <= imsize(1) );
        img_ind_ydir_neigh = sub2ind( imsize , img_ind_y( img_ind_ydir ) + 1 , img_ind_x( img_ind_ydir ) );

        ndir_enode1 = [ ndir_enode1 ; img_ind_ydir ];            
        ndir_enode2 = [ ndir_enode2 ; img_ind_ydir_neigh ];            
        ndir_ecap   = [ ndir_ecap   ; exp( -GradY(img_ind_ydir).^2 / VarNLink ) ];

        clear img_ind_ydir;
        clear img_ind_ydir_neigh;        
        
        % create sparse graph
        NeighGraph = sparse( ndir_enode1 , ...
                             ndir_enode2 , ... 
                             ndir_ecap , ...             
                             prod( imsize ) , prod( imsize ) ... 
                           );  
        
        % Display Spatial Cues
        figure, imshow( exp( -GradX.^2 / VarNLink ) , [] );
        title( 'Second-Order Clique Potentials -- XDir' );
        set( gcf , 'Name' , 'Second-Order Clique Potentials -- XDir' );         
        
        figure, imshow( exp( -GradY.^2 / VarNLink ) , [] );
        title( 'Second-Order Clique Potentials -- YDir' );
        set( gcf , 'Name' , 'Second-Order Clique Potentials -- YDir' );                 
        
        clear img_ind_*;
        
% Setup Unary-Clique Potentials

    cumnlinkcost = sum( ndir_ecap ) + 100;  
    
    % Bgnd
    FirstOrderCliquePotential_Bgnd = ( imInput(:) ~= 0 ) * cumnlinkcost;
    
    % Fgnd
    FirstOrderCliquePotential_Fgnd = ( imInput(:) == 0 ) * cumnlinkcost;
    
FirstOrderCliquePotential = [ FirstOrderCliquePotential_Bgnd , FirstOrderCliquePotential_Fgnd ];
    
% Solve Maxflow using KB maxflow wrapper
fprintf( 1 , '\n\n\tSolving the Maximum Flow Problem using KB maxflow wrapper ...\n\n'  );

    tic

    [ energy , seg_labels ] = maxflow( NeighGraph, sparse( cumnlinkcost - FirstOrderCliquePotential ) );

    energy

    imSegMask_kbmaxflow = reshape( seg_labels , imsize );
    imSegMask_kbmaxflow = imSegMask_kbmaxflow > 0;

    toc        
    
   figure, imshow( imSegMask_kbmaxflow , [] );
   title( 'Segmentation Result' );
   set( gcf , 'Name' , 'Segmentation Result' );     
   
% Check if test succeeded
if any( ( imSegMask_kbmaxflow(:) > 0 ) ~= ( imInput(:) > 0 ) ) 
    
    fprintf( 1 , '\n\nERROR: Test Failed ... \n\n' );
    
else
    
    fprintf( 1 , '\n\nTest Succeeded ... !!! \n\n' ); 
    
end

   