
function [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_strokes(im,displayrange)
% [varargout] = get_fgnd_bgnd_seeds_2d_strokes(im,range)
%
% This function allows the user to select get foreground and background
% seed points in a 2D slice using strokes
%

%%

if ~exist( 'displayrange' , 'var' )
    
    displayrange = [];
    
end

global data_get_fgnd_bgnd_seeds_2d_strokes;

hMainFigure = figure;
set( gcf , 'Name' , 'get_fgnd_bgnd_seeds_2d_strokes' );

% Create UI controls

    data_get_fgnd_bgnd_seeds_2d_strokes.ui.help_text = uicontrol(hMainFigure,'Style','text', ... 
                                                                'String','Provide Foreground and Background Seed Regions via 2D Strokes', ...
                                                                'BackgroundColor', [1 1 0.5], 'FontSize', 16, 'Units' , 'normalized' , 'Position',[0.15 0.92 0.4, 0.03]);  
    %[0.71 0.4 0.2, 0.1]);  
    
    % axis
    data_get_fgnd_bgnd_seeds_2d_strokes.ui.ah_img = axes( 'Position' , [ 0.001 , 0.2 , 0.7 , 0.7 ] , 'Visible' , 'off' );   
                
    % cursor point info controls
    data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_xloc = uicontrol(hMainFigure,'Style','edit',...
                    'String','X: INV',...
                    'Units' , 'normalized' , 'Position',[0.20 0.1 0.1 0.05]);                

    data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_yloc = uicontrol(hMainFigure,'Style','edit',...
                    'String','Y: INV',...
                    'Units' , 'normalized' , 'Position',[0.30 0.1 0.1 0.05]);     
                
    data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_Imval = uicontrol(hMainFigure,'Style','edit',...
                    'String','I: INV',...
                    'Units' , 'normalized' , 'Position',[0.40 0.1 0.1 0.05]);                                                   
    
    % seed selection mode controls       
    data_get_fgnd_bgnd_seeds_2d_strokes.ui.help_text = uicontrol(hMainFigure,'Style','text', ... 
                                                                'String','Seed Selection Mode: ', ...
                                                                'BackgroundColor', [1 1 1], 'FontSize', 14, 'Units' , 'normalized' , 'Position',[0.71 0.41 0.2, 0.03]);      
    data_get_fgnd_bgnd_seeds_2d_strokes.ui.bgh_mode = uibuttongroup('visible','on', 'Units' , 'normalized' ,'Position',[0.71 0.2 0.2 0.2]);
    data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_fgnd = uicontrol('Style','Radio','String','Foreground',...
                                 'Units' , 'normalized' ,'Position',[0.05 0.75 0.75 0.15],'parent',data_get_fgnd_bgnd_seeds_2d_strokes.ui.bgh_mode,'HandleVisibility','off');
    data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_bgnd = uicontrol('Style','Radio','String','Background',...
                                 'Units' , 'normalized' ,'Position',[0.05 0.50 0.75 0.15],'parent',data_get_fgnd_bgnd_seeds_2d_strokes.ui.bgh_mode,'HandleVisibility','off');            
    data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_none = uicontrol('Style','Radio','String','None',...
                                 'Units' , 'normalized' ,'Position',[0.05 0.25 0.75 0.15],'parent',data_get_fgnd_bgnd_seeds_2d_strokes.ui.bgh_mode,'HandleVisibility','off');    
    
    set( data_get_fgnd_bgnd_seeds_2d_strokes.ui.bgh_mode , 'SelectedObject' , data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_none );                             
    
% set callbacks
set( hMainFigure , 'WindowButtonDownFcn' , @FnMainFig_MouseButtonDownFunc );  
set( hMainFigure , 'WindowButtonMotionFcn' , @FnMainFig_MouseMotionFunc );  
set( hMainFigure , 'WindowButtonUpFcn' , @FnMainFig_MouseButtonUpFunc );  

% data_get_fgnd_bgnd_seeds_2d_strokes                         
data_get_fgnd_bgnd_seeds_2d_strokes.im = im;
data_get_fgnd_bgnd_seeds_2d_strokes.displayrange = displayrange;
data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points = {};
data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points = {};
data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging = false;

imsliceshow(data_get_fgnd_bgnd_seeds_2d_strokes);

% wait until the window is closed
errCatch = 0;
try
    waitfor( hMainFigure );
catch
    errCatch = 1;
end
    
if errCatch == 0 
    
    imsize = size( data_get_fgnd_bgnd_seeds_2d_strokes.im );

    fgnd_seed_points = [];
    
    if ~isempty( data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points )
    
        fgnd_seed_points = data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points;
        
        for cid = 1:numel( fgnd_seed_points )
            
            fgnd_seed_points{cid} = unique( round( fgnd_seed_points{cid} ) , 'rows' );
            
        end
        
    end
    
    bgnd_seed_points = [];
    
    if ~isempty( data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points )
    
        bgnd_seed_points = data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points;
        
        for cid = 1:numel( bgnd_seed_points )
            
            bgnd_seed_points{cid} = unique( round( bgnd_seed_points{cid} ) , 'rows' );
            
        end
        
    end
    
    clear data_get_fgnd_bgnd_seeds_2d_strokes;
    
else

    clear data_get_fgnd_bgnd_seeds_2d_strokes;
    error( 'Error: Unknown error occured while getting seed points from the user' );
    
end

    
    
%%
function imsliceshow(data_get_fgnd_bgnd_seeds_2d_strokes)

    imshow(data_get_fgnd_bgnd_seeds_2d_strokes.im,data_get_fgnd_bgnd_seeds_2d_strokes.displayrange);            

    hold on;

        if ~isempty( data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points )
            
            for cid = 1:numel( data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points )
                
                cur_fgnd_seed_points = data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points{cid};
                
                plot( cur_fgnd_seed_points( : , 1 ) , cur_fgnd_seed_points( : , 2 ) , '.g' );
                
            end

        end
        
        if ~isempty( data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points )
            
            for cid = 1:numel( data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points )
                
                cur_bgnd_seed_points = data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points{cid};
                
                plot( cur_bgnd_seed_points( : , 1 ) , cur_bgnd_seed_points( : , 2 ) , '.r' );
                
            end

        end
        
    hold off;
    
%%
function FnMainFig_MouseButtonDownFunc( hSrc , evnt )

    global data_get_fgnd_bgnd_seeds_2d_strokes;
    
    cp = get( gca , 'CurrentPoint' );
    
    if IsPointInsideImage( cp(1,1:2) , data_get_fgnd_bgnd_seeds_2d_strokes ) && strcmp( get(hSrc ,'SelectionType'),'normal' )       
        
        switch get( data_get_fgnd_bgnd_seeds_2d_strokes.ui.bgh_mode , 'SelectedObject' )
           
            case data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_fgnd
                
                data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points{end+1} = cp(1,1:2);

                data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging = true;
                
            case data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_bgnd
                
                data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points{end+1} = cp(1,1:2);
                
                data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging = true;                
        end       
        
    end
    
                   
    imsliceshow(data_get_fgnd_bgnd_seeds_2d_strokes);                
    

%% Update cursor point info -- xloc, yloc, int_val
function UpdateCursorPointInfo( data_get_fgnd_bgnd_seeds_2d_strokes )

    cp = get( gca , 'CurrentPoint' );       

    if IsPointInsideImage( cp(1,1:2) , data_get_fgnd_bgnd_seeds_2d_strokes )
        
        set(data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_xloc,'String' ,sprintf('X: %d / %d' , round( cp(1,1) ) , size( data_get_fgnd_bgnd_seeds_2d_strokes.im , 2 ) ));
        set(data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_yloc,'String' ,sprintf('Y: %d / %d' , round( cp(1,2) ) , size( data_get_fgnd_bgnd_seeds_2d_strokes.im , 1 ) ));        
        set(data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_Imval,'String',sprintf('I: %.1f' , data_get_fgnd_bgnd_seeds_2d_strokes.im( round( cp(1,2) ) , round( cp(1,1) ) ) ));                
        
    else
        
        set(data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_xloc,'String',sprintf('X: INV') );
        set(data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_yloc,'String',sprintf('Y: INV') );        
        set(data_get_fgnd_bgnd_seeds_2d_strokes.ui.eth_Imval,'String',sprintf('I: INV') );        
        
    end
    
    
%%    
function FnMainFig_MouseMotionFunc( hSrc , evnt )    
    
    global data_get_fgnd_bgnd_seeds_2d_strokes;
    
    cp = get( gca , 'CurrentPoint' );       
    
    if IsPointInsideImage( cp(1,1:2) , data_get_fgnd_bgnd_seeds_2d_strokes )
        
        set( hSrc ,'Pointer','crosshair');        
       
        if data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging 
        
            switch get( data_get_fgnd_bgnd_seeds_2d_strokes.ui.bgh_mode , 'SelectedObject' )

                case data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_fgnd

                    data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points{end} = [ data_get_fgnd_bgnd_seeds_2d_strokes.fgnd_seed_points{end} ; cp(1,1:2) ];

                case data_get_fgnd_bgnd_seeds_2d_strokes.ui_rbh_bgnd

                    data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points{end} = [ data_get_fgnd_bgnd_seeds_2d_strokes.bgnd_seed_points{end} ; cp(1,1:2) ];
                    
            end            
            
            imsliceshow(data_get_fgnd_bgnd_seeds_2d_strokes);
            
        end        
        
    else
        
        set( hSrc ,'Pointer','arrow');                
        
        if data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging
            
            data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging = false;
            
        end
        
    end        
                       
    UpdateCursorPointInfo( data_get_fgnd_bgnd_seeds_2d_strokes );

%%    
function FnMainFig_MouseButtonUpFunc( hSrc , evnt )

    global data_get_fgnd_bgnd_seeds_2d_strokes;
    
    cp = get( gca , 'CurrentPoint' );
    
    if IsPointInsideImage( cp(1,1:2) , data_get_fgnd_bgnd_seeds_2d_strokes ) && strcmp( get(hSrc ,'SelectionType'),'normal' )       
        
        if data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging
            
            data_get_fgnd_bgnd_seeds_2d_strokes.blnDragging = false;
            
        end
        
    end    
                       
  
%%    
function [ blnInside ] = IsPointInsideImage( cp , data_get_fgnd_bgnd_seeds_2d_strokes )

    imsize = size( data_get_fgnd_bgnd_seeds_2d_strokes.im );
    
    blnInside = all( cp <= imsize([2 1]) ) && all( cp >= [1 1] );
  
           