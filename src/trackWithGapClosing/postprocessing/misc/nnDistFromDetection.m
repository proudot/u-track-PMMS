function [nnDistMean,nnDistStd] = nnDistFromDetection(movieInfo)

nnDist = [];
for iFrame = 1 : length(movieInfo)
    
    if size(movieInfo(iFrame).xCoord,1) ~= 0
        
        xCoord1 = movieInfo(iFrame).xCoord(:,1);
        yCoord1 = movieInfo(iFrame).yCoord(:,1);
        
        featureDist = createDistanceMatrix([xCoord1 yCoord1],...
            [xCoord1 yCoord1]);
        
        featureDist = sort(featureDist,2);
        featureDist = featureDist(:,2);
        nnDist = [nnDist; featureDist];
        
    end
    
end

nnDistMean = mean(nnDist);
nnDistStd = std(nnDist);

%% ~~~ the end ~~~
