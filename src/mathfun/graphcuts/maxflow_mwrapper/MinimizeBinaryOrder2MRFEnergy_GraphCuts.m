function [ energy, labels ] = MinimizeBinaryOrder2MRFEnergy_GraphCuts( Order1CliquePotential , Order2Cliques , Order2CliquePotential )
% Minimized a given sub-modular (pairwise terms should be sub-modular) MRF-energy function of upto order-2 on binary variables using graph-cuts
% 
% Written By: Deepak Roy Chittajallu 
%
% This function is written based on the following paper of Kolmogorov et. al. :
% V. Kolmogorov and C. Rother, "Minimizing Nonsubmodular Functions with Graph Cuts-A Review,"
% IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 29, pp. 1274-1279, 2007.
% 
% Given a sub-modular energy function, this function does the following:
% 
% 1. Reduces/re-parameterizes the energy function to the normal form
% 2. Check for submodularity and check whether all (unary and pairwise) terms are non-negative
% 3. Constructs the graph for the energy function in the normal form
% 4. Minimizes the energy function using graph-cuts
%
% Input:  
% 
% Order1CliquePotential:
% 
%     It is a numNodes x numLabels matrix. 
%     
%     Order1CliquePotential( i , j ) measures the cost of assigning node i to label j
%
%     This function is currently limited to numLabels = 2 i.e. L = {0,1}  
%     
% Order2Cliques:
% 
%     It is a numEdges x 2 matrix. 
%     
%     Order2Cliques( i , : ) = ( from-node , to-node ) of edge i
%     
% Order2CliquePotential:
% 
%     It is a numLabels x numLabels x numEdges matrix. 
%     
%     Order2CliquePotential( L1 , L2 , i ) = Measures the cost of assigning the from-node and to-node of edge i 
%                                            to labels L1 and L2, respectively   
%     
% Output:
% 
%     energy: minimum MRF energy
%     
%     labels: minimum energy solution 
%

% validate inputs
if size( Order1CliquePotential , 2 ) ~= 2         
    
    error( 'ERROR: This function should be used to minimize submodular MRF-energies defined over binary variables only ...' );    
    
end

if size( Order2CliquePotential , 1 ) ~= 2 & size( Order2CliquePotential , 2 ) ~= 2 
    
    error( 'ERROR: This function should be used to minimize submodular MRF-energies defined over binary variables only ...' );    
    
end

if size( Order2Cliques , 1 ) ~= size( Order2CliquePotential , 3 )
    
    error( 'ERROR: Invalid Dimension of Order2Cliques and Order2CliquePotential ... They should have same number of edges' );    
    
end

Order1_pot = Order1CliquePotential;
Order2_pot = Order2CliquePotential;

numNodes = size( Order1_pot , 1 );
numEdges = size( Order2Cliques , 1 );

% Step-1: Reduce the unary and pairwise terms to normal form
E_const = 0;

    % reduce pairwise-terms
    for eind = 1:numEdges

        from_node_ind = Order2Cliques( eind , 1 );
        to_node_ind = Order2Cliques( eind , 2 );    
        pairwise_energies = Order2_pot( : , : , eind );

        for Lind = 1:2

            % normal-form check
            cur_min_val = min( pairwise_energies(:,Lind) );
            if cur_min_val ~= 0

                pairwise_energies(:,Lind) = pairwise_energies(:,Lind) - cur_min_val;
                Order1_pot( to_node_ind , Lind ) = Order1_pot( to_node_ind , Lind ) + cur_min_val;

            end

        end

        Order2_pot( : , : , eind ) = pairwise_energies;

        if any( pairwise_energies(:) < 0 )
            
            eind
            error( 'ERROR: pariwise term < 0 after reduction to normal form' );
            
        end
        
    end
    
    % reduce unary terms
    min_val_ldir = min( Order1_pot , [] , 2 );
    Order1_pot = Order1_pot - repmat( min_val_ldir , 1 , 2 );
    E_const = sum( min_val_ldir );

% Step-2: Check for submodularity and check whether all terms are non-negative

    % non-negative check
    if any( Order1_pot(:) < 0 ) | any( Order2_pot(:) < 0 ) | any( isnan( Order1_pot(:) ) ) | any( isnan( Order2_pot(:) ) )
    
        error( 'ERROR: Conversion to normal form failed -- after converting to normal form all terms (unary and pairwise) must be non-negative' );
        
    end
    
    % Fix Inf val (hard constraint) if present
    Infval = sum( Order1_pot( ~isinf(Order1_pot) ) ) + sum( Order2_pot( ~isinf(Order2_pot) ) ) + 100;
    
        % unary terms    
        if any( isinf( Order1_pot(:) ) )
            
            fprintf( 1 , '\nNote: Found Inf values in unary terms -- Fixing them to the highest possible value' );
            Order1_pot( isinf(Order1_pot) ) = Infval;
            
        end
            
        % pairwise terms
        if any( isinf( Order2_pot(:) ) )

            fprintf( 1 , '\nNote: Found Inf values in pairwise terms -- Fixing them to the highest possible value' );
            Order2_pot( isinf( Order2_pot ) ) = Infval;
        end
    
    % sub-modularity
    if any( squeeze( ~( Order2_pot( 1 , 1 , : ) + Order2_pot( 2 , 2 , : ) <= Order2_pot( 1 , 2 , : ) + Order2_pot( 2 , 1 , : ) ) ) )
    
        error( 'ERROR: The given energy function did not pass sub-modularity check' );
        
    end
    
% Step-3: Construct the graph for the energy function in the normal form
ndir_enode1 = Order2Cliques( : , 1 );
ndir_enode2 = Order2Cliques( : , 2 );

NeighGraph = sparse( [ ndir_enode1 ; ndir_enode2 ] , ... % from-node
                     [ ndir_enode2 ; ndir_enode1 ] , ... % to-node
                     [ squeeze( Order2_pot( 1 , 2 , : ) ) ; squeeze( Order2_pot( 2 , 1 , : ) ) ] , ... % edge capacities
                     numNodes , numNodes );

% Step-4: Minimize the energy function using graph-cuts
[ min_energy , labels ] = maxflow( NeighGraph, sparse( fliplr( Order1_pot ) ) );

E_const 

%energy = E_const + min_energy;
energy = min_energy;

end