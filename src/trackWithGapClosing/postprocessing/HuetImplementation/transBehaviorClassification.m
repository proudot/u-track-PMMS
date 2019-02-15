function trajClassification = transBehaviorClassification(classParameters,...
                             trajectoryMatrix, tWinSize, nDiff, nDev, nMSD);
% TRANSBEHAVIORCLASSIFICATION function which classifieds the
%                             transient behavior of trajectories into 4
%                             classes and plots them colorcoded:
%                             stalled(cyan), constrained(red), directed 
%                             (green), brownian(blue) and unclassified(black) 
%
% SYNOPSIS trajClassification = transBehaviorClassification(classParameters,...
%                               trajectoryMatrix, tWinSize, nDiff, nDev, nMSD)
%
% INPUT  classParameters  = 3-by-2 matrix with 
%                             rows: 1. parameters for stalled behavior
%                                   2. parameters for constrained behavior
%                                   3. parameters for directed behavior
%                             columns: 1. thresholds for parameters
%                                      2. crossing times for classification
%
%        trajectoryMatrix = n-by-8 matrix with
%                             n:       trajectories
%                             columns: x, y, intensity, amplitude,
%                                      dx, dx, dintensity, damplitude
%                             
%        tWinSize         = 3-by-2 matrix with 
%                             rows: 1. for diffusion
%                                   2. for constrained
%                                      min win size should be 1.5 times
%                                      nDev
%                                   3. for directed
%                             columns : 1. minimum Size of the rolling time
%                                          window
%                                       2. maximum size of the rolling time
%                                          window       
%                           Guidelines for choosing the time window size:
%                         - the minimum size of the time window should be
%                           smaller then the minimum duration of a certain
%                           behavior. This can only be rigorously verified
%                           a posteriori to the analysis.
%                         - the mean duration of a trajectory could be used
%                           as the maximum window size. It can be reduced
%                           by a priori estimation with a test pool of 
%                           trajectories to save computational time
%
%        nDiff (optional) = number of timpoints of MSD which are used for
%                           the linear regression fit to determine the 
%                           diffussion coeffitions. (default = 5)
%
%        nDev (optional)  = number of timepoints of MSD compared to the
%                           linear regression fit to determine the
%                           deviation between them. (default = 50)
%
%        nMSD (optional)  = number of timepoints used for the calculation
%                           of the MSD (min = 4, max = tWinSizeMin-2,
%                           default = 20) 
%
% OUTPUT  trajClassification  = n-by-18-by-l matrix with
%             n rows:  timepoints of trajectories
%             columns: 1. minimum diffusion coeffitient of all different 
%                         time windows width for each time point
%                      2. p-Value, which reflects the probability that
%                         sigma Hat Square is fischer distributed
%                         with a mean of 1.
%                      3. minimum deviation of the MSD from linearity (fit 
%                         over first nDiff points of MSD) of all different 
%                         time windows width for each time point
%                      4. maximum asymmetry of all different time windows 
%                         width for each time point 
%                      5. classification of transient behavior
%                                 0 = stalled
%                                 1 = brownian 
%                                 2 = constrained
%                                 3 = directed
%                      6. constrainded classified == 2
%                      7. directed classified == 3
%                      8. window size used for computing the min diffusion
%                         coeffitient in column 1
%                      9. window size used for computing the min deviation
%                         of the MSD from linearity in column 3
%                     10. window size used for computing the max asymmetry
%                         in column 3
%                  11-15. classification extended over the window size
%                         (column 8-10) used for computation of parameters
%                         (column 1-4)
%                               11. stalled classified == 0
%                               12. brownian classified == 1
%                               13. constrained classified == 2
%                               14. directed classified == 3
%                               15. final classification
%                                       0 = stalled
%                                       1 = brownian 
%                                       2 = constrained
%                                       3 = directed  
%                     16. actual number of data points used for the
%                         calculation of the min diffusion coeffitient
%                         (column 1)
%                     17. actual number of data points used for the
%                         calculation of the min deviation of the MSD from
%                         linearity (column 3)
%                     18. actual number of data points used for the
%                         calculation of the max asymmetry (column 4)
%             l: trajectories
%
% CREATED  gp 4/11/07

%-------------------------------------------
% initialize
%-------------------------------------------

trajClassification = [];

% ----------------------------------------------
% checking the input variables
% ----------------------------------------------

% checks if the number of time points to calculate the MSD (nMSD) is given
% by the user, if not it is set to default(20)
if (nargin <= 5);
    nMSD = 20;
end

% checks if the number of time points to calculate the deviation (nDev) is
% given by the user, if not it is set to default(5)
if (nargin <= 4);
    nDev = 50;
end

% checks if the number of time points to calculate the diffusion (nDiff) is
% given by the user, if not it is set to default(5)
if (nargin <= 3);
    nDiff = 5;
end

% calles the function transBehaviorAnalysis which analysis the transient 
% behavior of trajectories by calculating the diffusion coeffitient, the
% deviation of the MSD and the assymetry with a rolling time window of
% varying size for the different parameters.
[trajClassification] = transBehaviorAnalysis(trajectoryMatrix, tWinSize,...
 nDiff, nDev, nMSD);

% deternines the size of the matrix trajClassification
[rowSize, columnSize, thirdDimension] = size(trajClassification);

% determines the rows where all parameters are available for the decission
% of the type of movement
allParameters = find(sum(isnan(trajClassification(:,1:4))')' == 0);



%--------------------------------------------------------------
% classification of behavior charakteristics of each timepoint
%--------------------------------------------------------------

% loops through the different trajectories 
for n = 1:thirdDimension
    
    % calculation of devRel which determines the relative deviation
    % corresponding to the min window size which was used to calculate the
    % minimum deviation of each time point 
    trajClassification(:,3,n) = trajClassification(:,3,n) - ((-0.0097*(log(trajClassification(:,9,n)).^2))...
        + 0.1413*(log(trajClassification(:,9,n))) -1.1085);
 

    % detection of directed periods: all positions which show an asymmetry
    % above the treshhold for directed periods are defined as directed (3)
    directed = find(trajClassification(:,4,n) >= classParameters(3,1));
    trajClassification(directed,7,n) = 3;
    trajClassification(directed,5,n) = 3;
    
    
    % detection of constrained periods: all positions which show a deviation
    % below the treshhold for constrained periods are defined as constrained (2)
    constrained = find(trajClassification(:,3,n) <= classParameters(2,1));
    trajClassification(constrained,6,n) = 2;
    trajClassification(constrained,5,n) = 2;
    
    % detection of stalled periods: all positions which show a diffusion
    % coeffitient below the treshhold for stalled periods are defined as
    % stalled (0)
    stalled = find(trajClassification(:,1,n) <= classParameters(1,1));
    trajClassification(stalled,5,n) = 0;
     
    
    
            
%     % detection of brownian periods: all positions which are not classified so
%     % far will be determined as browian (1)
%     isnan(trajClassification(min(allParameters):max(allParameters),5,n));
%     for i = min(allParameters):max(allParameters)
%         if isnan(trajClassification(i,5,n));
%            trajClassification(i,5,n) = 1;
%         end
%     end
%     
    %-------------------------------------------------------------
    % detects directed periods
    %-------------------------------------------------------------
    
    % all directed periods are set to 1 and all others to 0 in the
    % vector findDir
    findDir = trajClassification(:,5,n) == 3;
 
    % sets the first and the last timepoint to 0 to find beginning and ends
    % of classifications when the reach the ends
    findDir(1) = 0;
    findDir(end) = 0;
   
    % substracts n+1 - n which will result in a vector containing 1 for the
    % start and -1 for the end of a directed period
    findDir = diff(findDir);

    % executes only if there is a directed area
    if any(findDir) == 1

        % finds the indices of the beginnings of the directed periods
        dirStarts = (find(findDir == 1))+1;

        % corrects the start of the classified period if its at the very first
        % point of the trajectory by substracting 1
        if dirStarts(1) == 2
            dirStarts(1) = 1;
        end
        
        % finds the indices of the ends of the directed periods
        dirEnds = find(findDir == -1);

        % corrects the end of the classified period if its at the very last
        % point of the trajectory by adding 1
        if dirEnds(end) == size(trajClassification,1)-1
            dirEnds(end) = dirEnds(end) + 1;
        end
       
        % determines the length of the directed periods and compares them
        % with the set crossing time for directed movements from the
        % classParameters
        ind = find(dirEnds-dirStarts >= classParameters(3,2));
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % uncomment the whole paragraph to get no expansion of the
%         % directed movement
%
%         for i = 1:size(ind,1)
% 
%            % setting the classifaction to directed in clolumn 15 without
%            % expanding the directed movement over the timewindow used
%            % for determination of the max asymmertry
%            trajClassification((dirStarts(ind(i))+1):(dirEnds(ind(i))),15,n) = 3;
% 
%         end % of for i = 1:size(ind,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % uncomment the whole paragraph to get expansion of the
        % directed movement
        
        for i = 1:size(ind,1)
            tWinSizeStart = trajClassification(dirStarts(ind(i))+1,10);
            tWinSizeEnd = trajClassification(dirEnds(ind(i)),10);
            movementWinStart =  dirStarts(ind(i))+ 1 - ((tWinSizeStart - 1)/2);
            movementWinEnd =  dirEnds(ind(i)) + ((tWinSizeEnd - 1)/2);
            trajClassification(movementWinStart:movementWinEnd,14,n) = 3;
            trajClassification(movementWinStart:movementWinEnd,15,n) = 3;  
        end % of for i = 1:size(ind,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    end % of if any(findDir) == 1
        
    %-------------------------------------------------------------
    % detects constrained periods
    %-------------------------------------------------------------
    
    % all constrained periods are set to 1 and all others to 0 in the
    % vector findcon
    findCon = trajClassification(:,5,n) == 2;
    
    % sets the first and the last timepoint to 0 to find beginning and ends
    % of classifications when the reach the ends
    findCon(1) = 0;
    findCon(end) = 0;
   
    % substracts n+1 - n which will result in a vector containing 1 for the
    % start and -1 for the end of a constrained period
    findCon = diff(findCon);
   
    % executes only if there is a constrained area
    if any(findCon) == 1
       
        % finds the indices of the beginnings of the constrained periods
        conStarts = (find(findCon == 1))+1;
       
        % corrects the start of the classified period if its at the very first
        % point of the trajectory by substracting 1
        if conStarts(1) == 2
            conStarts(1) = 1;
        end
       
        % finds the indices of the ends of the constrained periods
        conEnds = find(findCon == -1);
     
        % corrects the end of the classified period if its at the very last
        % point of the trajectory by adding 1
        if conEnds(end) == (size(trajClassification,1))-1
            conEnds(end) = (conEnds(end)) + 1;
        end

        % determines the length of the constrained periods and compares them
        % with the set crossing time for constrained movements from the
        % classParameters
        ind = find(conEnds-conStarts >= classParameters(2,2));
    
        for i = 1:size(ind,1)
            tWinSizeStart = trajClassification(conStarts(ind(i))+1,9);
            tWinSizeEnd = trajClassification(conEnds(ind(i)),9);
            movementWinStart =  conStarts(ind(i))+ 1 - ((tWinSizeStart - 1)/2);
            movementWinEnd =  conEnds(ind(i)) + ((tWinSizeEnd - 1)/2);
            trajClassification(movementWinStart:movementWinEnd,13,n) = 2;
            trajClassification(movementWinStart:movementWinEnd,15,n) = 2;
        end % of for i = 1:size(ind,1)
        
    end % of if any(findCon) == 1

% %......................................................................... 
% %finding overlaps between directed and constrained motion and dividing the
% %overlaps into 2/3 constrained and 1/3 directed 
% %.........................................................................
%     if any((trajClassification(:,13,n) == 2) & (trajClassification(:,14,n) == 3));
%         %find overlap between direted and constrained classification
%         conDirOverLap = (trajClassification(:,13,n) == 2) +...
%             (trajClassification(:,14,n) == 3);
%         conDirOverLap(conDirOverLap==0)=1;
%         conDirOverLap(1) = 1;
%         conDirOverLap = diff(conDirOverLap);
%         conDirOverLapStarts = find(conDirOverLap == 1)+1;
%         conDirOverLapStarts(conDirOverLapStarts==2)=1;
%         conDirOverLapEnds = find(conDirOverLap == -1);
% 
%         for i=1:length(conDirOverLapEnds)
% 
%             % lenght of overlap
%             overLapLength = conDirOverLapEnds(i)-conDirOverLapStarts(i);
%             
%             % if classification was constrained before overlap
%             if conDirOverLapStarts(i) ~= 1
%                 if trajClassification(conDirOverLapStarts(i)-1,13,n) == 2
%                     % 2/3 of overlaplength are classified constrained
%                     trajClassification(conDirOverLapStarts(i):(floor(2/3*overLapLength)),15,n) = 2;
%                 % if classification was directed before overlap    
%                 elseif trajClassification(conDirOverLapStarts(i)-1,14,n) == 3
%                     % 1/3 of overlaplength are classified directed
%                     trajClassification(conDirOverLapStarts(i):conDirOverLapStarts(i)+(floor(1/3*overLapLength)),15,n) = 3;
%                 end
%             end
%             
%              % if classification was constrained after overlap
%             if trajClassification(conDirOverLapEnds(i)+1,13,n) == 2
%                 % 2/3 of overlaplength are classified constrained
%                 trajClassification((conDirOverLapEnds(i)-(floor(2/3*overLapLength))):conDirOverLapEnds(i),15,n) = 2;
%             % if classification was directed after overlap    
%             elseif trajClassification(conDirOverLapEnds(i)+1,14,n) == 3
%                 % 1/3 of overlaplength are classified directed
%                 trajClassification((conDirOverLapEnds(i)-(floor(1/3*overLapLength))):conDirOverLapEnds(i),15,n) = 3;
%             end
%             
%         end
% %........................................................................

    %-------------------------------------------------------------
    % detects stalled periods
    %-------------------------------------------------------------
    
    % all stalled periods are set to 1 and all others to 0 in the
    % vector findStalled
    findStalled = trajClassification(:,5,n) == 0;
    
    
        
    % sets the first and the last timepoint to 0 to find beginning and ends
    % of classifications when the reach the ends
    findStalled(1) = 0;
    findStalled(end) = 0;
    
    % substracts n+1 - n which will result in a vector containing 1 for the
    % start and -1 for the end of a stalled period
    findStalled = diff(findStalled);
    
    % executes only if there is a constrained area
    if any(findStalled) == 1
        
        % finds the indices of the beginnings of the stalled periods
        stalledStarts = find(findStalled == 1);

        % finds the indices of the ends of the stalled periods
        stalledEnds = find(findStalled == -1);

        % determines the length of the stalled periods and compares them
        % with the set crossing time for stalled movements from the
        % classParameters
        ind = find(stalledEnds-stalledStarts >= classParameters(1,2));

        for i = 1:size(ind,1)
            tWinSizeStart = trajClassification(stalledStarts(ind(i))+1,8);
            tWinSizeEnd = trajClassification(stalledEnds(ind(i)),8);
            movementWinStart =  stalledStarts(ind(i))+ 1 - ((tWinSizeStart - 1)/2);
            movementWinEnd =  stalledEnds(ind(i)) + ((tWinSizeEnd - 1)/2);
            trajClassification(movementWinStart:movementWinEnd,11,n) = 0;
            trajClassification(movementWinStart:movementWinEnd,15,n) = 0;
        end % of i = 1:size(ind,1)
    end % of if any(findStalled) == 1
    %-------------------------------------------------------------
    % detects brownian movement
    %-------------------------------------------------------------
    
    % detection of brownian periods: all positions which are not classified so
    % far will be determined as browian (1)
    isnan(trajClassification(min(allParameters):max(allParameters),15,n));
    for i = min(allParameters):max(allParameters)
        if isnan(trajClassification(i,15,n));
           trajClassification(i,15,n) = 1;
        end
    end

   %--------------------------------------------------------------
   % plotting of different transient behaviors
   %--------------------------------------------------------------
   
   % opens a new figure for every trajectory
   figure;
   
   % keeps points in graph by adding new ones
   hold on;
   
   % plot the whole trajectory in black
   plot(trajectoryMatrix(n,1:8:end), trajectoryMatrix(n,2:8:end),'k');
              
   %-------------------------------------------------------------
   % plotting of constrained segments
   %-------------------------------------------------------------

   % all constrained periods are set to 1 and all others to 0 in the
   % vector findcon
   findCon = (trajClassification(:,15,n) == 2);
   
   % sets the first and the last timepoint to 0 to find beginning and ends
   % of classifications when the reach the ends
   findCon(1) = 0;
   findCon(end) = 0;
   
   % substracts n+1 - n which will result in a vector containing 1 for the
   % start and -1 for the end of a constrained period
   findCon = diff(findCon);
   
   % executes only if there is a constrained area
   if any(findCon) == 1
       
       % finds the indices of the beginnings of the constrained periods
        conStarts = (find(findCon == 1))+1;
       
       % corrects the start of the classified period if its at the very first
       % point of the trajectory by substracting 1
       if conStarts(1) == 2
           conStarts(1) = 1;
       end
       
       % finds the indices of the ends of the constrained periods
       conEnds = find(findCon == -1);
     
       % corrects the end of the classified period if its at the very last
       % point of the trajectory by adding 1
       if conEnds(end) == (size(trajClassification,1))-1
           conEnds(end) = (conEnds(end)) + 1;
       end
       
       % loops through the different segments of constrained movement
       for i = 1:size(conStarts)

           % plots the constrained parts of the trajectory in magenta
           plot(trajectoryMatrix(n,(conStarts(i)+1)*8-7:8:conEnds(i)...
               *8-7),trajectoryMatrix(n,(conStarts(i)+1)*8-6:8:conEnds(i)...
               *8-6),'r'); 

       end % of constrained plotting
   end % of if any(findCon) == 1

   %-------------------------------------------------------------
   % plotting of stalled segments
   %-------------------------------------------------------------

   % all stalled periods are set to 1 and all others to 0 in the
   % vector findStalled
   findStalled = (trajClassification(:,15,n) == 0);
   
   % sets the first and the last timepoint to 0 to find beginning and ends
   % of classifications when the reach the ends
   findStalled(1) = 0;
   findStalled(end) = 0;
   
   % substracts n+1 - n which will result in a vector containing 1 for the
   % start and -1 for the end of a stalled period
   findStalled = diff(findStalled);
   
   % executes only if there is a stalled area
   if any(findStalled) == 1
       
        % finds the indices of the beginnings of the stalled periods
       stalledStarts = (find(findStalled == 1))+1;

       % corrects the start of the classified period if its at the very first
       % point of the trajectory by substracting 1
       if stalledStarts(1) == 2
           stalledStarts(1) = 1;
       end
       
       % finds the indices of the ends of the stalled periods
       stalledEnds = find(findStalled == -1);
     
       % corrects the end of the classified period if its at the very last
       % point of the trajectory by adding 1
       if stalledEnds(end) == (size(trajClassification,1))-1
           stalledEnds(end) = (stalledEnds(end)) + 1;
       end

       % loops through the different segments of stalled movement
       for i = 1:size(stalledStarts)

           % plots the stalled parts of the trajectory in green
           plot(trajectoryMatrix(n,(stalledStarts(i)+1)*8-7:8:stalledEnds(i)...
                *8-7),trajectoryMatrix(n,(stalledStarts(i)+1)*8-6:8:stalledEnds(i)...
                *8-6),'c'); 

       end % of stalled plotting
   end % of if any(findStalled) == 1 

   %-------------------------------------------------------------
   % plotting of directed segments
   %-------------------------------------------------------------

   % all directed periods are set to 1 and all others to 0 in the
   % vector findcon
   findDir = (trajClassification(:,15,n) == 3);
   
   % sets the first and the last timepoint to 0 to find beginning and ends
   % of classifications when the reach the ends
   findDir(1) = 0;
   findDir(end) = 0;
   
   % substracts n+1 - n which will result in a vector containing 1 for the
   % start and -1 for the end of a directed period
   findDir = diff(findDir);

   % executes only if there is a directed area
   if any(findDir) == 1

       % finds the indices of the beginnings of the directed periods
       dirStarts = (find(findDir == 1))+1;

       % corrects the start of the classified period if its at the very first
       % point of the trajectory by substracting 1
       if dirStarts(1) == 2
           dirStarts(1) = 1;
       end
       
       % finds the indices of the ends of the directed periods
       dirEnds = find(findDir == -1);

       % corrects the end of the classified period if its at the very last
       % point of the trajectory by adding 1
       if dirEnds(end) == size(trajClassification,1)-1
           dirEnds(end) = dirEnds(end) + 1;
       end

       % loops through the different segments of directed movement
       for i = 1:size(dirStarts)

           % plots the directed parts of the trajectory in green
           plot(trajectoryMatrix(n,(dirStarts(i))*8-7:8:dirEnds(i)...
               *8-7),trajectoryMatrix(n,(dirStarts(i))*8-6:8:dirEnds(i)...
               *8-6),'g'); 

        end % of directed plotting
   end % of if any(findDir) == 1

   %-------------------------------------------------------------
   % plotting of brownian segments
   %-------------------------------------------------------------

   % all brownian periods are set to 1 and all others to 0 in the
   % vector findcon
   findBrown = (trajClassification(:,15,n) == 1);
   
   % sets the first and the last timepoint to 0 to find beginning and ends
   % of classifications when the reach the ends
   findBrown(1) = 0;
   findBrown(end) = 0;
   
   % substracts n+1 - n which will result in a vector containing 1 for the
   % start and -1 for the end of a brownian period
   findBrown = diff(findBrown);
   
   % executes only if there is a brownian area
   if any(findBrown) == 1

       % finds the indices of the beginnings of the brownian periods
       brownStarts = find(findBrown == 1);
       
       % corrects the start of the classified period if its at the very first
       % point of the trajectory by substracting 1
       if brownStarts(1) == 2
           brownStarts(1) = 1;
       end
       
       % finds the indices of the ends of the brownian periods
       brownEnds = find(findBrown == -1);
       
       % corrects the end of the classified period if its at the very last
       % point of the trajectory by adding 1
       if brownEnds(end) == size(trajClassification,1)
           brownEnds(end) = brownEnds(end)+1; 
       end
       
       % loops through the different segments of brownian movement
       for i = 1:size(brownStarts)

           % plots the brownian parts of the trajectory in cyan
           plot(trajectoryMatrix(n,(brownStarts(i)+1)*8-7:8:brownEnds(i)...
               *8-7),trajectoryMatrix(n,(brownStarts(i)+1)*8-6:8:brownEnds(i)*8-6),'b'); 

       end % of brownian plotting
   end % of if any(findBrown) == 1 
      
end










