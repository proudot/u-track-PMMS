This folder contains the matlab code that implements simple interactive segmentation using graph-cuts

The main function/entry point is segObject2D_maxflow() which can be found in the file segObject2D_maxflow.m

Run the script test_segObject2D_maxflow.m to see a demo on a synthetically generated image.

If u want to apply the function on some other image, you will have to add the subdirectories "maxflow_mwrapper" and "maxflow_mwrapper\maxflow-v3.0" to the matlab path. 
Also you will have to run "maxflow_mwrapper/make_maxflow.m" to generate the mex files of the maxflow algorithm. You can do this with the following commands:

addpath('maxflow_mwrapper', 'maxflow_mwrapper/maxflow-v3.0');
cd maxflow_mwrapper
make_maxflow
cd ..

 
Function Usage and Description:

[ imSegMask_kbmaxflow ] = segObject2D_maxflow( imInput , DisplayRange )

	segObject2D_maxflow allows the user to interactively segment an image into foreground and background regions

	This function implements the traditional interactive image segmentation using graph-cuts (See Ref[1] below). 

	As soon as the function is called, a GUI window will popup asking the user to interactively provide the foreground and background seed regions via 2D 
	strokes. 
	
	The strokes provided using this GUI are used to build appearance models for the foreground and the background object. 
	Both foreground and background appearance is modeled using a mixture of gaussians. Also, one gaussian/normal pdf is constructed for each stroke provided 
	by the user.

	After obtaining the seed regions, a 2-label MRF problem is setup and solved to obtain the segmentation.

	Sites of the MRF are pixels in the image
	Labels - {foregound, background}
	Cliques - We only consider first- and second-order cliques. The second-order cliques are defined based on a 4-neighborhood graph. 
			  So, there are two types of order-2 cliques: (i) between neighboring pixels in x-direction, and (ii) between neighboring pixels in y-direction

	The first-order clique potentials corresponding to foreground and background are determined using the mixture of guassian models built from user-define 
	seed regions. 

	V_i( L_i = 1 ) = -log( pdf( I(i) , mu_fgnd ; var_fgnd ) );
	V_i( L_i = 0 ) = -log( pdf( I(i) , mu_bgnd ; var_bgnd ) );

	The second-order clique potentials are defined as a generalized potts interaction model based on a function that is inversely proportional to the 
	intensity gradient between the two pixels involved in the clique.

	V_ij( L_i , L_j ) = 0 if L_i = L_j
					  = exp( -( I(i) - I(j) )^2 / (2 * sigma^2) ) if L_i ~= L_j 

	As might be evident, the above defined second-order clique is submodular and hence can be minimized using graph-cuts.
	We use the implementation of Boykov and Komogrov to solve the maxflow problem (http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html).

Input arguments:

--  imInput: 

    A 2D matrix containing the grayscale image that needs to be segmented.
    
    If imInput is not provided then segmentation is performed on a synthetically generated image. This gives you a demo of what the function does and how it 
	works.
    
-- DisplayRange (optional): 

   An intensity range for better visualization of the input image. 
   Default value [min(imInput(:) max(imInput(:))]   
   
Output arguments:

--  imSegMask_kbmaxflow

    A binary segmentation mask

Example usages:

-- Demo Mode

    imMask = segObject2D_maxflow();    
    
-- Example-1

    I = imread('coins.png');
    seg_mask = segObject2D_maxflow( I , [0 255] );

References:

[1] Y. Boykov and G. Funka-Lea. Graph cuts and efficient ND image segmentation. International Journal of Computer Vision, 70(2):109–131, 2006.
[2] V. Kolmogorov and R. Zabin. What energy functions can be minimized via graph cuts? IEEE Trans. on Pattern Analysis and Machine Intelligence, 
	26(2):147–159, 2004.
[3] Y. Boykov, V. Kolmogorov. An experimental comparison of min-cut/max-flow algorithms for energy minimization in vision. IEEE Trans Pattern Analysis and 
	Machine Intelligence, 26:1124–1137, 2004

