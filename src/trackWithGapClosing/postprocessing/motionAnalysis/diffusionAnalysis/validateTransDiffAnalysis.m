function [matrixCompPerPart,matrixCompPerPoint,fracPointsCorrectPerPart,...
    aveStartEndOffset,startOffsetDistr,endOffsetDistr] = ...
    validateTransDiffAnalysis(classGT,classTest)

%initialize the results matrices
matrixCompPerPart = zeros(4);
matrixCompPerPoint = zeros(4);
fracPointsCorrectPerPart = zeros(4,3);
aveStartEndOffset = zeros(4);
startOffsetDistr = struct('confined',[],'brownian',[],'directed',[]);
endOffsetDistr = startOffsetDistr;

%get number of tracks
numTracks = length(classGT);

%go over all tracks
for iTrack = 1 : numTracks

    %get number of segments in current track
    numSegments = length(classGT(iTrack).segmentClass);
    
    %go over all segments
    for iSegment = 1 : numSegments
    
        %get segment classification from ground truth and from transient
        %diffusion analysis
        segmentClassGT = classGT(iTrack).segmentClass(...
            iSegment).momentScalingSpectrum(:,1:3);
        segmentClassTest = classTest(iTrack).segmentClass(...
            iSegment).momentScalingSpectrum(:,1:3);

        %replace NaN, indicating unclassified, with 0, to allow its
        %counting later
        segmentClassGT(isnan(segmentClassGT(:,3)),3) = 0;
        segmentClassTest(isnan(segmentClassTest(:,3)),3) = 0;
        
        %get number of classification subparts from ground truth and from
        %transient diffusion analysis
        numClassPartsGT = size(segmentClassGT,1);
        numClassPartsTest = size(segmentClassTest,1);
        
        %get number of time points in track segment
        nTP = segmentClassGT(end,2);
        
        %make an array of point classifications from ground truth and from
        %transient diffusion analysis
        pointClassGT = NaN(nTP,1);
        for iPart = 1 : numClassPartsGT
            pointClassGT(segmentClassGT(iPart,1):segmentClassGT(iPart,2)) = ...
                segmentClassGT(iPart,3);
        end
        pointClassTest = NaN(nTP,1);
        for iPart = 1 : numClassPartsTest
            pointClassTest(segmentClassTest(iPart,1):segmentClassTest(iPart,2)) = ...
                segmentClassTest(iPart,3);
        end

        %COMPARISON 1:
        %go over track segment subparts from the ground truth and compare
        %to classification from transient diffusion analysis
        for iPart = 1 : numClassPartsGT
            
            %get this subpart's GT information
            partStartGT = segmentClassGT(iPart,1);
            partEndGT = segmentClassGT(iPart,2);
            partClassGT = segmentClassGT(iPart,3);
            
            %get the correponsing array of point classifications from
            %transient diffusion analysis
            correspondPointClassArrayTest = pointClassTest(partStartGT:partEndGT);
            
            %get the number of points in each classification
            %0 = unclassified, 1 = confined, 2 = Brownian, 3 = directed
            numPointsEachClassTest = hist(correspondPointClassArrayTest,(0:3));
            
            %find the classification with the largest number of points
            %this is considered to be the classification of this part from
            %the transient diffusion analysis, corresponding to the
            %ground truth classification
            partClassTest = find( numPointsEachClassTest == ...
                max(numPointsEachClassTest) ) - 1;
            
            %take care of the case when two classses have the same number
            %of points by chance ...
            if length(partClassTest) > 1 && any(partClassTest == partClassGT)
                partClassTest = partClassGT;
            else
                partClassTest = partClassTest(1);
            end
                
            
            %add 1 to the entry representing this combination of ground
            %truth and transient diffusion analysis classifications
            matrixCompPerPart(partClassGT+1,partClassTest+1) = ...
                matrixCompPerPart(partClassGT+1,partClassTest+1) + 1;
            
            %if the two classifications agree and this is neither the first
            %part nor the last part ...
            if partClassTest == partClassGT && iPart > 1 && iPart < numClassPartsGT
                
                %SUB-COMPARISON 1:
                %store the length of the subpart from the GT and the 
                %fraction of points with the correct classification from
                %the transient diffusion analysis
                
                %get part length and fraction statistics so far
                prevAveLength = fracPointsCorrectPerPart(partClassGT+1,1);
                prevFracCorrect = fracPointsCorrectPerPart(partClassGT+1,2);
                numPrevObs = fracPointsCorrectPerPart(partClassGT+1,3);
                
                %get current part length and fraction
                currentLength = partEndGT - partStartGT + 1;
                currentFracCorrect = max(numPointsEachClassTest) / currentLength;

                %combine the two
                fracPointsCorrectPerPart(partClassGT+1,:) = [...
                    mean([prevAveLength*ones(1,numPrevObs) currentLength]) ...
                    mean([prevFracCorrect*ones(1,numPrevObs) currentFracCorrect]) ...
                    numPrevObs+1];

                %SUB-COMPARISON 2:
                %store the offset between the real start and end of the
                %part and the start and end found by the transient
                %diffusion analysis

                %find the offest of the start of this part according to the
                %transient diffusion analysis
                partStartTestOffset = find(correspondPointClassArrayTest==partClassGT,1,'first') - 1;

                %if the start is at the first GT time point of this part,
                %go back and check whether the part starts earlier
                %according to the transient diffusion analysis
                if partStartTestOffset == 0
                    pointClassBeforeInverted = pointClassTest(partStartGT-1:-1:1);
                    partStartTestOffset = -find(pointClassBeforeInverted~=partClassGT,1,'first') + 1;
                    if isempty(partStartTestOffset)
                        partStartTestOffset = 1 - partStartGT;
                    end
                end

                %find the offest of the end of this part according to the
                %transient diffusion analysis
                partEndTestOffset = -find(correspondPointClassArrayTest(end:-1:1)==partClassGT,1,'first') + 1;

                %if the end is at the last GT time point of this part, go
                %forward and check whether the part ends earlier according
                %to the transient diffusion analysis
                if partEndTestOffset == 0
                    pointClassAfterNotInverted = pointClassTest(partEndGT+1:end);
                    partEndTestOffset = find(pointClassAfterNotInverted~=partClassGT,1,'first') - 1;
                    if isempty(partEndTestOffset)
                        partEndTestOffset = length(pointClassAfterNotInverted);
                    end
                end

                %get offset statistics so far
                prevStartOffsetNP = aveStartEndOffset(partClassGT+1,1);
                prevStartOffsetPP = aveStartEndOffset(partClassGT+1,2);
                prevEndOffsetNP   = aveStartEndOffset(partClassGT+1,3);
                prevEndOffsetPP   = aveStartEndOffset(partClassGT+1,4);

                %combine the previous offset statistics with the current
                %offset statistics
                aveStartEndOffset(partClassGT+1,:) = [...
                    mean([prevStartOffsetNP*ones(1,numPrevObs) partStartTestOffset]) ...
                    mean([prevStartOffsetPP*ones(1,numPrevObs) abs(partStartTestOffset)]) ...
                    mean([prevEndOffsetNP*ones(1,numPrevObs) partEndTestOffset]) ...
                    mean([prevEndOffsetPP*ones(1,numPrevObs) abs(partEndTestOffset)])];
                
                %store alse the offset distribution
                switch partClassGT
                    case 1
                        startOffsetDistr.confined = [startOffsetDistr.confined; ...
                            partStartTestOffset];
                        endOffsetDistr.confined = [endOffsetDistr.confined; ...
                            partEndTestOffset];
                    case 2
                        startOffsetDistr.brownian = [startOffsetDistr.brownian; ...
                            partStartTestOffset];
                        endOffsetDistr.brownian = [endOffsetDistr.brownian; ...
                            partEndTestOffset];
                    case 3
                        startOffsetDistr.directed = [startOffsetDistr.directed; ...
                            partStartTestOffset];
                        endOffsetDistr.directed = [endOffsetDistr.directed; ...
                            partEndTestOffset];
                end

            end %(if partClassTest == partClassGT)

        end %(for iPart = 1 : numClassPartsGT)

        %COMPARISON 2:
        %go over each point and compare its ground truth classification to
        %its classification from transient diffusion analysis
        for iClassGT = 0 : 3
            for iClassTest = 0 : 3
                matrixCompPerPoint(iClassGT+1,iClassTest+1) = ...
                    matrixCompPerPoint(iClassGT+1,iClassTest+1) + ...
                    length(find( pointClassGT == iClassGT & ...
                    pointClassTest == iClassTest ));
            end
        end

    end %(for iSegment = 1 : numSegments)

end %(for iTrack = 1 : numTracks)

%for matrixCompPerPart,  sum up all the numbers, then get
%fractions relative to ground truth
sumGT = sum(matrixCompPerPart,2);
sumTest = sum(matrixCompPerPart,1);
matrixCompPerPart = matrixCompPerPart ./ repmat(sumGT,1,4);
matrixCompPerPart = [matrixCompPerPart sumGT];
matrixCompPerPart(isnan(matrixCompPerPart)) = 0;
matrixCompPerPart = [matrixCompPerPart; [sumTest NaN]];

%do the same for matrixCompPerPoint
sumGT = sum(matrixCompPerPoint,2);
sumTest = sum(matrixCompPerPoint,1);
matrixCompPerPoint = matrixCompPerPoint ./ repmat(sumGT,1,4);
matrixCompPerPoint = [matrixCompPerPoint sumGT];
matrixCompPerPoint(isnan(matrixCompPerPoint)) = 0;
matrixCompPerPoint = [matrixCompPerPoint; [sumTest NaN]];

%remove the last column (with number of observations) from
%fracPointsCorrectPerPart
fracPointsCorrectPerPart = fracPointsCorrectPerPart(:,1:2);

