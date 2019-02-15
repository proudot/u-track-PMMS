
function [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_points(im,displayrange)
% [varargout] = get_fgnd_bgnd_seeds_2d_points(im,range)
%
% This function allows the user to select get foreground and background
% seed points in a 2D slice using points
%

%%

if ~exist( 'displayrange' , 'var' )
    
    displayrange = [];
    
end

global data_get_fgnd_bgnd_seeds_2d_points;

hMainFigure = figure;

% Create UI controls

    % axis
    data_get_fgnd_bgnd_seeds_2d_points.ui.ah_img = axes( 'Position' , [ 0.001 , 0.2 , 0.7 , 0.7 ] , 'Visible' , 'off' );
   
    % cursor point info controls
    data_get_fgnd_bgnd_seeds_2d_points.ui.eth_xloc = uicontrol(hMainFigure,'Style','edit',...
                    'String','X: INV',...
                    'Units' , 'normalized' , 'Position',[0.20 0.05 0.1 0.05]);                

    data_get_fgnd_bgnd_seeds_2d_points.ui.eth_yloc = uicontrol(hMainFigure,'Style','edit',...
                    'String','Y: INV',...
                    'Units' , 'normalized' , 'Position',[0.30 0.05 0.1 0.05]);     
                
    data_get_fgnd_bgnd_seeds_2d_points.ui.eth_Imval = uicontrol(hMainFigure,'Style','edit',...
                    'String','I: INV',...
                    'Units' , 'normalized' , 'Position',[0.40 0.05 0.1 0.05]);                                                
                
    % seed selection mode controls
    data_get_fgnd_bgnd_seeds_2d_points.ui.bgh_mode = uibuttongroup('visible','on', 'Units' , 'normalized' ,'Position',[0.71 0.2 0.2 0.2]);
    data_get_fgnd_bgnd_seeds_2d_points.ui_rbh_fgnd = uicontrol('Style','Radio','String','Foreground',...
                                 'Units' , 'normalized' ,'Position',[0.05 0.75 0.75 0.15],'parent',data_get_fgnd_bgnd_seeds_2d_points.ui.bgh_mode,'HandleVisibility','off');
    data_get_fgnd_bgnd_seeds_2d_points.ui_rbh_bgnd = uicontrol('Style','Radio','String','Background',...
                                 'Units' , 'normalized' ,'Position',[0.05 0.50 0.75 0.15],'parent',data_get_fgnd_bgnd_seeds_2d_points.ui.bgh_mode,'HandleVisibility','off');            
    data_get_fgnd_bgnd_seeds_2d_points.ui_rbh_none = uicontrol('Style','Radio','String','None',...
                                 'Units' , 'normalized' ,'Position',[0.05 0.25 0.75 0.15],'parent',data_get_fgnd_bgnd_seeds_2d_points.ui.bgh_mode,'HandleVisibility','off');    
    
    set( data_get_fgnd_bgnd_seeds_2d_points.ui.bgh_mode , 'SelectedObject' , data_get_fgnd_bgnd_seeds_2d_points.ui_rbh_none );                             
    
% set callbacks
set( hMainFigure , 'WindowButtonDownFcn' , @FnMainFig_MouseButtonDownFunc );  
set( hMainFigure , 'WindowButtonMotionFcn' , @FnMainFig_MouseMotionFunc );  

% data_get_fgnd_bgnd_seeds_2d_points                         
data_get_fgnd_bgnd_seeds_2d_points.im = im;
data_get_fgnd_bgnd_seeds_2d_points.displayrange = displayrange;
data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points = [];
data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points = [];

imsliceshow(data_get_fgnd_bgnd_seeds_2d_points);

% wait until the window is closed
errCatch = 0;
try
    waitfor( hMainFigure );
catch
    errCatch = 1;
end
    
if errCatch == 0 
    
    imsize = size( data_get_fgnd_bgnd_seeds_2d_points.im );

    fgnd_seed_points = [];
    
    if ~isempty( data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points )
    
        fgnd_seed_points = unique( round( data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points ) , 'rows' );
        
    end
    
    bgnd_seed_points = [];
    
    if ~isempty( data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points )
    
        bgnd_seed_points = unique( round( data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points ) , 'rows' );
        
    end
    
    clear data_get_fgnd_bgnd_seeds_2d_points;
    
else
    
    clear data_get_fgnd_bgnd_seeds_2d_points;
    error( 'Error: Unknown error occured while getting seed points from the user' );
    
end
    
%%
function imsliceshow(data_get_fgnd_bgnd_seeds_2d_points)

    imshow(data_get_fgnd_bgnd_seeds_2d_points.im,data_get_fgnd_bgnd_seeds_2d_points.displayrange);            

    hold on;

        if ~isempty( data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points )
            
            plot( data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points( : , 1 ) , data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points( : , 2 ) , 'g+' );

        end
        
        if ~isempty( data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points )
            
            plot( data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points( : , 1 ) , data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points( : , 2 ) , 'r+' );

        end
        
    hold off;

%%
function FnMainFig_MouseButtonDownFunc( hSrc , evnt )

    global data_get_fgnd_bgnd_seeds_2d_points;
    
    cp = get( gca , 'CurrentPoint' );
    
    if IsPointInsideImage( cp(1,1:2) , data_get_fgnd_bgnd_seeds_2d_points ) && strcmp( get(hSrc ,'SelectionType'),'normal' )       
        
        switch get( data_get_fgnd_bgnd_seeds_2d_points.ui.bgh_mode , 'SelectedObject' )
           
            case data_get_fgnd_bgnd_seeds_2d_points.ui_rbh_fgnd
                
                data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points = [ data_get_fgnd_bgnd_seeds_2d_points.fgnd_seed_points ; cp(1,1:2) ];

            case data_get_fgnd_bgnd_seeds_2d_points.ui_rbh_bgnd
                
                data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points = [ data_get_fgnd_bgnd_seeds_2d_points.bgnd_seed_points ; cp(1,1:2) ];                
        end

    end
    
    imsliceshow(data_get_fgnd_bgnd_seeds_2d_points);
    

%% Update cursor point info -- xloc, yloc, int_val
function UpdateCursorPointInfo( data_get_fgnd_bgnd_seeds_2d_points )
    
    cp = get( gca , 'CurrentPoint' );       

    if IsPointInsideImage( cp(1,1:2) , data_get_fgnd_bgnd_seeds_2d_points )
        
        set(data_get_fgnd_bgnd_seeds_2d_points.ui.eth_xloc,'String' ,sprintf('X: %d / %d' , round( cp(1,1) ) , size( data_get_fgnd_bgnd_seeds_2d_points.im , 2 ) ));
        set(data_get_fgnd_bgnd_seeds_2d_points.ui.eth_yloc,'String' ,sprintf('Y: %d / %d' , round( cp(1,2) ) , size( data_get_fgnd_bgnd_seeds_2d_points.im , 1 ) ));        
        set(data_get_fgnd_bgnd_seeds_2d_points.ui.eth_Imval,'String',sprintf('I: %.1f' , data_get_fgnd_bgnd_seeds_2d_points.im( round( cp(1,2) ) , round( cp(1,1) ) ) ));                
        
    else
        
        set(data_get_fgnd_bgnd_seeds_2d_points.ui.eth_xloc,'String',sprintf('X: INV') );
        set(data_get_fgnd_bgnd_seeds_2d_points.ui.eth_yloc,'String',sprintf('Y: INV') );        
        set(data_get_fgnd_bgnd_seeds_2d_points.ui.eth_Imval,'String',sprintf('I: INV') );        
        
    end

%%    
function FnMainFig_MouseMotionFunc( hSrc , evnt )    
       
    global data_get_fgnd_bgnd_seeds_2d_points;

    cp = get( gca , 'CurrentPoint' );       
    
    if IsPointInsideImage( cp(1,1:2) , data_get_fgnd_bgnd_seeds_2d_points )
        
        set( hSrc ,'Pointer','crosshair');        
        
    else
        
        set( hSrc ,'Pointer','arrow');        

    end    
    
    
    imsliceshow(data_get_fgnd_bgnd_seeds_2d_points);    
    UpdateCursorPointInfo( data_get_fgnd_bgnd_seeds_2d_points );

%%    
function [ blnInside ] = IsPointInsideImage( cp , data_get_fgnd_bgnd_seeds_2d_points )

    imsize = size( data_get_fgnd_bgnd_seeds_2d_points.im );
    
    blnInside = all( cp <= imsize([2 1]) ) && all( cp >= [1 1] );
  

            