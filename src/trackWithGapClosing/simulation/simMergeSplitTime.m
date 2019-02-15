function [mergeTimeDistr,splitTimeDistr] = simMergeSplitTime(probDim,...
    diffCoef,totalTime,timeStep,numRepeat,resLimit)

%initialize output vectors
mergeTimeDistr = NaN(numRepeat,1);
splitTimeDistr = NaN(numRepeat,1);

%repeat numRepeat times
for i=1:numRepeat

    %generate the position time series of particle 1
    pos1 = brownianMotion(probDim,diffCoef,totalTime+1,timeStep);
    
    %generate the position time series of particle 2, such that the initial
    %distance between particles 1 and 2 is the resolution limit
    pos2 = resLimit/sqrt(probDim) + brownianMotion(probDim,diffCoef,totalTime+1,timeStep);
    
    %average over intervals of 1 time unit
    pos1Ave = mean(reshape(pos1(1:end-1,1),1/timeStep,[]))';
    pos2Ave = mean(reshape(pos2(1:end-1,1),1/timeStep,[]))';
    if probDim > 1
        y1Ave = mean(reshape(pos1(1:end-1,2),1/timeStep,[]))';
        y2Ave = mean(reshape(pos2(1:end-1,2),1/timeStep,[]))';
        pos1Ave = [pos1Ave y1Ave];
        pos2Ave = [pos2Ave y2Ave];
    end
    if probDim > 2
        z1Ave = mean(reshape(pos1(1:end-1,3),1/timeStep,[]))';
        z2Ave = mean(reshape(pos2(1:end-1,3),1/timeStep,[]))';
        pos1Ave = [pos1Ave z1Ave];
        pos2Ave = [pos2Ave z2Ave];
    end

    %calculate the distance between the two particles
    dist = sqrt(sum((pos1Ave - pos2Ave).^2,2));
    dist = dist(2:end);

    if dist(1) < resLimit %if the particles merged ...

        %find the earliest time where the distance between the two
        %particles becomes larger than the resolution limit
        timeM = find(dist>resLimit,1,'first');

        %if they do split, enter the time of merging
        if ~isempty(timeM)
            mergeTimeDistr(i) = timeM - 1;
        else
            mergeTimeDistr(i) = totalTime;
        end

    elseif dist(1) > resLimit %if the particles split

        %find the earliest time where the distance between the two
        %particles becomes smaller than the resolution limit
        timeS = find(dist<resLimit,1,'first');

        %if they do merge, enter the time of splitting
        if ~isempty(timeS)
            splitTimeDistr(i) = timeS - 1;
        else
            splitTimeDistr(i) = totalTime;
        end

    end %(if dist(2) < resLimit ... elseif ...)

end

