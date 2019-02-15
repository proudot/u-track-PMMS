function y = nancumsumwithzero(x)

%find non-NaN values
indxNonNaN = find(~isnan(x));

%extract consecutive intervals of non-NaNs
indxDiff = diff(indxNonNaN);
indxJump = find(indxDiff > 1);
intervalStart = [indxNonNaN(1); indxNonNaN(indxJump+1)];
intervalEnd = [indxNonNaN(indxJump); indxNonNaN(end)];

%get number of intervals
numInterval = length(intervalStart);

%go over each interval and calculate the cumulative sum
y = NaN(size(x));
for iInt = 1 : numInterval
    y(intervalStart(iInt):intervalEnd(iInt)) = cumsum(x(intervalStart(iInt):intervalEnd(iInt)));
end

%add zeros so that diff(y) = x
y(intervalStart(2:end)-1) = 0;
y(2:end+1) = y;
y(1) = 0;
