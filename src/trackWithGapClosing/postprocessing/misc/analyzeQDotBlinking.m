
%define imaging conditions to be compared
exposure = [8 16 30];
sensitivity = [150 200 240];
laserPower = [1 2 3 4];

for i = 1 : length(exposure)
    for j = 1 : length(sensitivity)
        for k = 1 : length(laserPower)

            %get directory name
            directoryName = ['/mnt/sickkids/Yoav/2009_05_06_monodisperse_QDots_blinking_calibration/' ...
                num2str(exposure(i)) 'ms/' num2str(sensitivity(j)) 'sens/cell_0' num2str(laserPower(k)) '/'];
            analysisDir = [directoryName 'analysis/'];

            %load tracks
            cd(analysisDir)
            load tracksTest3_1.mat
            
            %convert tracks into matrix format
            tracks = convStruct2MatNoMS(tracksFinal);
            
            %find tracks 4 frames or longer
            criteria.lifeTime.min = 4;
            indx4 = chooseTracks(tracks,criteria);
            tracks4 = tracks(indx4,:);
            numTracks4 = length(indx4);
            
            %calculate number of features per frame from tracks and
            %tracks4, hence calculate fraction of false positives
            xCoord = tracks(:,1:8:end);
            xCoord4 = tracks4(:,1:8:end);
            numFeatures = zeros(100,1);
            numFeatures4 = numFeatures;
            fracFalsePos = numFeatures;
            for l = 1 : 100
                numFeatures(l) = length(find(~isnan(xCoord(:,l))));
                numFeatures4(l) = length(find(~isnan(xCoord4(:,l))));
                fracFalsePos(l) = 1 - (numFeatures4(l) / numFeatures(l));
            end
            
            %calculate the average of the above quantities
            aveNumFeat = mean(numFeatures);
            aveNumFeat4 = mean(numFeatures4);
            aveFracFalsePos = mean(fracFalsePos);
            
            %get track start, end and life times
            trackSEL4 = getTrackSEL(tracks4);
            
            %get track gaps
            trackGaps4 = findTrackGaps(tracks4);
            
            %calculate number of gaps, average gap length and maximum gap
            %length
            numGaps4 = size(trackGaps4,1);
            aveGapLength4 = mean(trackGaps4(:,4));
            maxGapLength4 = max(trackGaps4(:,4));
            
            %calculate average number of gaps per feature
            aveNumGapsPerFeat4 = numGaps4 / aveNumFeat4;
            
            %calculate overall intensity mean and std
            intensityVec = tracks4(:,4:8:end);
            intensityVec = intensityVec(intensityVec~=0&~isnan(intensityVec));
            intensityMean4 = mean(intensityVec);
            intensityStd4 = std(intensityVec);
            
            %calculate the average localization precision in x
            xStd4 = tracks4(:,5:8:end);
            xStd4 = xStd4(~isnan(xStd4));
            imagTmp = imag(xStd4);
            xStd4 = real(xStd4(imagTmp==0));
            xStd4 = xStd4(xStd4<1000);
            aveStdX = mean(xStd4);

            %calculate the average localization precision in y
            yStd4 = tracks4(:,6:8:end);
            yStd4 = yStd4(~isnan(yStd4));
            imagTmp = imag(yStd4);
            yStd4 = real(yStd4(imagTmp==0));
            yStd4 = yStd4(yStd4<1000);
            aveStdY = mean(yStd4);


            %store all variables in the output variable
            analysis(k,j,i) = struct('tracks',tracks,'tracks4',tracks4,...
                'numFeatures',numFeatures,'numFeatures4',numFeatures4,...
                'fracFalsePos',fracFalsePos,'aveNumFeat',aveNumFeat,...
                'aveNumFeat4',aveNumFeat4,'aveFracFalsePos',aveFracFalsePos,...
                'trackGaps4',trackGaps4,'numGaps4',numGaps4,...
                'aveGapLength4',aveGapLength4,'maxGapLength4',maxGapLength4,...
                'aveNumGapsPerFeat4',aveNumGapsPerFeat4,'intensityMean4',...
                intensityMean4,'intensityStd4',intensityStd4,'aveStdX',...
                aveStdX,'aveStdY',aveStdY);

        end
    end
end

%extract some information out of the analysis structure and save results
directoryName = '/mnt/sickkids/Yoav/2009_05_06_monodisperse_QDots_blinking_calibration/';
cd(directoryName);
save('QdotAnalysisRes','analysis','exposure','sensitivity','laserPower');
clear all
load QdotAnalysisRes
for i = 1 : length(exposure)
    for j = 1 : length(sensitivity)
        for k = 1 : length(laserPower)
            intensityMean(k,j,i) = analysis(k,j,i).intensityMean4;
            intensityStd(k,j,i) = analysis(k,j,i).intensityStd4;
            fracFalsePos(k,j,i) = analysis(k,j,i).aveFracFalsePos;
            posStdXY(k,j,i) = mean([analysis(k,j,i).aveStdX analysis(k,j,i).aveStdY]);
            numGapsPerFeat(k,j,i) = analysis(k,j,i).aveNumGapsPerFeat4;
            aveGapLength(k,j,i) = analysis(k,j,i).aveGapLength4;
            maxGapLength(k,j,i) = analysis(k,j,i).maxGapLength4;
            
        end
    end
end
clear i j k
save('QdotAnalysisRes')


%% Old stuff

%             %for every track, calculate its intensity mean and standard
%             %deviation, as well as its gap characteristics
%             indIntMean = zeros(numTracks4,1);
%             indIntStd = indIntMean;
%             indFracGaps = indIntMean;
%             indGapLength = indIntMean;
%             for l = 1 : numTracks4
%                 intensityVec = tracks4(l,4:8:end);
%                 indGapsVec = trackGaps4(trackGaps4(:,1)==l,4);
%                 indIntMean(l) = nanmean(intensityVec);
%                 indIntStd(l) = nanstd(intensityVec);
%                 indFracGaps(l) = length(indGapsVec)/trackSEL4(l,3);
%                 indGapLength(l) = mean(indGapsVec);
%             end

