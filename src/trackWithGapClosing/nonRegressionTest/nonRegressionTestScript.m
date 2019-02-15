% A QD script to compare the performances of multiple tracking branches (tested with 2 ...)
% Warning: this script uses  simulated data as a benchmark. Please synchronise simulation function in each branch before testing. 


%% Branches setup
% Each branch must be cloned in a separated directory
branchNb=2;
branches(branchNb) = struct();

branches(1).name='iu-track';
branches(1).path='/home2/proudot/repo/iu-track-branch/';

branches(2).name='master';
branches(2).path='/home2/proudot/repo/danuser-utsw/';

arrayfun(@(x) (rmpath(genpath(x.path))),branches);

%% simulation parameter
% 1 for classic simulation (brownian or confined ordirected)
% 2 for heterogeneous
simulation_type=1;

% Varying param for simulation density in spots/1 microm^2 (20 px window)
speed_transition_range=1:10;

renderSimulation=false;

%% tracking methods
% Tracking method(s) compared
method_nb=1;
compute_U_track=1;
compute_KF_iter_mostRecent=0;


%% Memory allocation for results for each simulation and each tracking method
for i=1:length(branches)
    branches(i).percentageGoodLink=zeros(method_nb,size(speed_transition_range,2),4);
    branches(i).percentageBadLink=zeros(method_nb,size(speed_transition_range,2),4);
    branches(i).percentageMissingLink=zeros(method_nb,size(speed_transition_range,2),4);
    branches(i).noise_var_MSE=zeros(method_nb,size(speed_transition_range,2));
    branches(i).cpuTime=zeros(method_nb,size(speed_transition_range,2));
    branches(i).schemesStats=zeros(2,method_nb,size(speed_transition_range,2));
end

%% For each parameter in the input range, simulated and track for each branch.
simuParamIdx=1;
for speed_transition=speed_transition_range
    addpath(genpath(branches(1).path));

    %% TracksSim parametrization and computation
    detectionParam.bitDepth = 16; %Camera bit depth
    imSize=[100 100];
    density=5;
    numP=density*(imSize(1)*imSize(2))/400
    %numP=50;
    numF=100;
    sigmaLft=20;
    meanLft=50;
    x =1:numF;
    lftDist=exp(-(x-meanLft).^2/(2*sigmaLft^2))/sqrt(2*pi*sigmaLft^2);
    intVec=[500 30];
    
    switch simulation_type
      case 1
        motionParam.diffCoef2D=[0.1 4];
        motionParam.confRad2D=[0.5 2];
        motionParam.speed1D=[speed_transition-0.5 speed_transition+0.5];
        motionParam.probVelSwitch=0.8;
        motionParam.probMS=[];
        motionParam.fracType=[0.5 0.4 0.1];
        [simMPM,tracksSim] = simulateRandomPlusLinearMotion(imSize,numP,lftDist,numF,intVec,motionParam);
        trackDiffCoef=ones(1,length(tracksSim));  % FIX : get from simu
      case 2 
        motionParam.diffCoef2D=[0.1 2];
        motionParam.confRad2D=[0.5 2];
        motionParam.speed1D=[speed_transition-0.5 speed_transition+0.5];
        motionParam.probVelSwitch=0.5;
        motionParam.probTypeSwitch=0.1;
        motionParam.probMS=[];
        motionParam.fracType=[0 0 0 0 0 0 1];        
        [simMPM,tracksSim,trackDiffCoef] = simulateRandomPlusLinearMotionHetero(imSize,numP,lftDist,numF,intVec,motionParam);
    end 
    
    %% data translation
    % detection GT to movie info data structure
    % MovieInfo data structure allows to assign an implicite ID to each spot
    % We have to store this tracks in the ID form to link kalman filter info to real GT tracks property.
    movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numF,1);
    for i = 1: numF
        xCoord=simMPM(simMPM(:,3*i-2)~=1.,3*i-2);
        movieInfo(i).xCoord=[ xCoord 0.05*ones(size(xCoord))];
        yCoord=simMPM(simMPM(:,3*i-1)~=1.,3*i-1);
        movieInfo(i).yCoord=[ yCoord 0.05*ones(size(yCoord))];    
        amp=simMPM(simMPM(:,3*i)~=0.,3*i)/(2^detectionParam.bitDepth -1);
        movieInfo(i).amp=[ amp 0.0003*ones(size(amp))];        
    end

    % transcripting movieInfo ID on tracksSIm
    cumID = cumsum(simMPM(:,1:3:end)~=1,1);
    for track_idx = 1:size(tracksSim,1)
        tracksSim(track_idx).tracksFeatIndxCG=cumID(track_idx,simMPM(track_idx,1:3:end)~=1);
    end 

    %% simulated stack production (opt)
    if(renderSimulation)
        bgav=100;                               
        bgnoise=30;
        sigma=1.;
        rad=4;
        saveVar=1;
        saveFolder=['../../simulation/density_variation/density_' num2str(density) '/'];
        mkdir(saveFolder);
        imsize=imSize;
        revSimMPM=simMPM;                       % solving strange issue with makeAiryImageFromMPM
        revSimMPM(:,1:3:end)=simMPM(:,2:3:end);
        revSimMPM(:,2:3:end)=simMPM(:,1:3:end);
        makeAiryImageFromMPM(revSimMPM,bgav,bgnoise,sigma,imsize,rad,saveVar,saveFolder);
    end


    saveResults.filename = 'detectionTestnoGMM.mat'; %name of file where input and output are saved
    probDim=2;
    verbose =1;

    saveFolder='/tmp/';
    %% Tracker parametrization and computation
    movieParam.imageDir = saveFolder;
    movieParam.filenameBase = 'simulImage'; %image file name base
    movieParam.firstImageNum = 1; %number of first image in movie
    movieParam.lastImageNum =numF; %number of last image in movie
    movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).

    [gapCloseParam,costMatrices,kalmanFunctions,probDim,verbose]=control_KF_tracker_param();

    method_idx=1;              

    %% U track
    if compute_U_track
        
        schemeName='U_track'; %directory where to save input and output
        
        %function name
        costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';
        kalmanFunctions.reserveMem  = 'kalmanResMemLM';
        kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
        kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
        kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
        
        %recording results
       
        for i=1:length(branches)
            addpath(genpath(branches(i).path));
            [branches(i).percentageGoodLink(method_idx,simuParamIdx,:), ...
             branches(i).percentageBadLink(method_idx,simuParamIdx,:), ...
             branches(i).percentageMissingLink(method_idx, simuParamIdx,:), ... 
             branches(i).noise_var_MSE(method_idx,simuParamIdx),...
             branches(i).schemesStats(:,method_idx, simuParamIdx),...
             branches(i).cpuTime(method_idx,simuParamIdx)] = ...
                track_stat(movieInfo,tracksSim, ...
                           trackDiffCoef,costMatrices,kalmanFunctions,gapCloseParam,saveFolder,movieParam, ...
                           schemeName);
            rmpath(genpath(branches(i).path)); 

        end
        method_idx=method_idx+1;
        method_idx
        
    end

    %% KF iter
    if compute_KF_iter_mostRecent
        
        schemeName='KF_iter_mostRecent';
        
        %function name
        costMatrices(1).funcName = 'costMatRandDirSwitchIndepKalmanLinkIter';
        kalmanFunctions.reserveMem  = 'indeKalmanResMemLMOnlineVar';
        kalmanFunctions.initialize  = 'indeKalmanInitLinearMotionOnlineVar';
        kalmanFunctions.calcGain    = 'indeKalmanGainLinearMotionRebDirectNoise';
        kalmanFunctions.timeReverse = 'indeKalmanReverseLinearMotion';
        
        for i=1:length(branches)
            addpath(genpath(branches(i).path));
            [branches(i).percentageGoodLink(method_idx,simuParamIdx,:), ...
             branches(i).percentageBadLink(method_idx,simuParamIdx,:), ...
             branches(i).percentageMissingLink(method_idx, simuParamIdx,:), ... 
             branches(i).noise_var_MSE(method_idx,simuParamIdx),...
             branches(i).schemesStats(:,method_idx, simuParamIdx),...
             branches(i).cpuTime(method_idx,simuParamIdx)] = ...
                track_stat(movieInfo,tracksSim, ...
                           trackDiffCoef,costMatrices,kalmanFunctions,gapCloseParam,saveFolder,movieParam, ...
                           schemeName);
            rmpath(genpath(branches(i).path)); 

        end
       
        method_idx=method_idx+1;
        method_idx
    end 
    
    simuParamIdx=simuParamIdx+1;
end

save('branches-performance.mat','branches')


%% plotting results
figure(1);
cc=hsv(branchNb);

subplot(2,2,1);
for i=1:length(branches)
plot(speed_transition_range,branches(i).percentageGoodLink(:,:,1)','+-','color',cc(i,:));
hold on
end
hold off
axis([min(speed_transition_range),max(speed_transition_range),0,100])
title('U-track linking percentage')
legend(branches(:).name)
xlabel('transition speed')
ylabel('correct link percentage')

subplot(2,2,2);
for i=1:length(branches)
    plot(speed_transition_range,branches(i).percentageBadLink(:,:,1)','+-','color',cc(i,:));
    hold on
end
hold off
axis([min(speed_transition_range),max(speed_transition_range),0,100])
title('U-track wrong linking percentage')
legend(branches(:).name)
xlabel('transition speed')
ylabel('wrong link percentage')

subplot(2,2,4);
for i=1:length(branches)
    plot(speed_transition_range,branches(i).percentageBadLink(:,:,1)','+-','color',cc(i,:));
hold on
end
hold off
axis([min(speed_transition_range),max(speed_transition_range),0,100])
title('U-track linking percentage')
legend(branches(:).name)
xlabel('Motion type switching probability')
ylabel('Good link type 4 ???')
title('linking percentage on direct motion') 

subplot(2,2,3);
for i=1:length(branches)
    plot(speed_transition_range,sqrt(branches(i).noise_var_MSE)','+-','color',cc(i,:));
hold on
end
hold off
title('U-track process noise RMSE')
legend(branches(:).name)
xlabel('transition speed')
ylabel('Process noise error (px)')

%% Results summary
figure(3)
for i=1:length(branches)
    plot(speed_transition_range,branches(i).percentageGoodLink(:,:,1)','+-','color',cc(i,:));
end
axis([min(speed_transition_range),max(speed_transition_range),0,100])
legend(branches(:).name)
xlabel('transition speed')
ylabel('correct link percentage')


%% FP to
subplot(1,2,1);
for i=1:length(branches)
    plot(speed_transition_range,branches(i).percentageGoodLink(:,:,1)','+-','color',cc(i,:));
hold on
end
hold off
axis([min(speed_transition_range),max(speed_transition_range),0,100])
legend(branches(:).name)
xlabel('transition speed')
ylabel('correct link percentage')

subplot(1,2,2);
for i=1:length(branches)
    plot(speed_transition_range,branches(i).percentageBadLink(:,:,1)','+-','color',cc(i,:));
hold on
end
hold off
axis([min(speed_transition_range),max(speed_transition_range),0,100])
legend(branches(:).name)
xlabel('transition speed')
ylabel('wrong link percentage')


%% check diff
resDiff(branchNb-1)=struct();
for i=1:(branchNb-1)
    resDiff(i).percentageGoodLink=branches(1).percentageGoodLink - branches(i+1).percentageGoodLink;
    resDiff(i).percentageBadLink=branches(1).percentageBadLink - branches(i+1).percentageBadLink;
end

if ( all(arrayfun(@(x) (~any(x.percentageGoodLink(:))),resDiff)) && all(arrayfun(@(x) (~any(x.percentageBadLink(:))),resDiff)) ) 
    disp('No regression detected');
else 
    disp('ERROR: please check resDiff');
end 




