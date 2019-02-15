%% Build pdf from seed points
function [ model ] = BuildModelPdf( seed_points , volInput , SEED_NEIGH )

    model = [];
        
    if ~isempty( seed_points )
        
        volSize = size( volInput );
        
        if iscell( seed_points )
            
            model.im_seed_mask = zeros( volSize ) > 0;

            for cid = 1:numel( seed_points )

                seed_points{cid} = unique( round( seed_points{cid} ) , 'rows' );        

                cur_seed_points = seed_points{cid};

                cur_im_seed_mask = zeros( volSize(1:2) );

                if size(cur_seed_points, 2) == 2
                    cur_sliceno = 1;
                else
                    cur_sliceno = cur_seed_points(1,3);
                end
                cur_dicom_slice = volInput( : , : , cur_sliceno );           

                cur_ptind = sub2ind( volSize(1:2) , cur_seed_points(:,2) , cur_seed_points(:,1) );

                for i = 1:numel( cur_ptind );

                    cx = cur_seed_points(i,1);
                    cy = cur_seed_points(i,2);

                    nx = max( 1 , cx - SEED_NEIGH ) : min( volSize(2) , cx + SEED_NEIGH );        
                    ny = max( 1 , cy - SEED_NEIGH ) : min( volSize(1) , cy + SEED_NEIGH );

                    cur_im_seed_mask( ny , nx ) = 1;

                end

                model.mean(cid,:) = mean( cur_dicom_slice(cur_im_seed_mask > 0) );
                model.covariance(:,:,cid) = cov( cur_dicom_slice(cur_im_seed_mask > 0) );              
                model.tissue_gmobj{cid} = gmdistribution( model.mean(cid,:) , model.covariance(:,:,cid) );

                model.im_seed_mask( : , : , cur_sliceno ) = model.im_seed_mask( : , : , cur_sliceno ) | cur_im_seed_mask;            

            end

            model.seed_points = seed_points;
            model.gmobj = gmdistribution( model.mean , model.covariance );
            
        else
           
            model_pdf.im_seed_mask = zeros( volSize ) > 0;

            cur_seed_points = unique( round( seed_points) , 'rows' );

            model.im_seed_mask = zeros( volSize );
            
            for i = 1:size( cur_seed_points , 1 )                       

                if size(cur_seed_points, 2) == 2
                    cur_sliceno = 1;
                else
                    cur_sliceno = cur_seed_points(i,3);
                end
                
                cur_dicom_slice = volInput( : , : , cur_sliceno );
                cur_im_seed_mask = zeros( volSize(1:2) );

                cx = cur_seed_points(i,1);
                cy = cur_seed_points(i,2);

                nx = max( 1 , cx - SEED_NEIGH ) : min( volSize(2) , cx + SEED_NEIGH );        
                ny = max( 1 , cy - SEED_NEIGH ) : min( volSize(1) , cy + SEED_NEIGH );

                cur_im_seed_mask( ny , nx ) = 1;

                model.mean(i,:) = mean( cur_dicom_slice(cur_im_seed_mask > 0) );
                model.covariance(:,:,i) = cov( cur_dicom_slice(cur_im_seed_mask > 0) );              
                model.sed_gmobj = gmdistribution( model.mean(i,:) , model.covariance(:,:,i) );

                model.im_seed_mask( : , : , cur_sliceno ) = model.im_seed_mask( : , : , cur_sliceno ) | cur_im_seed_mask;            

            end

            model.seed_points = cur_seed_points;
            model.gmobj = gmdistribution( model.mean , model.covariance );
            
        end
        
    end
