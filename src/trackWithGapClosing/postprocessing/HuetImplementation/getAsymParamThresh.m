
numSamples = 2000;
numTimePoints = 101;

traj1 = NaN*ones(numTimePoints,2,numSamples);
traj0p1 = NaN*ones(numTimePoints,2,numSamples);
asym1 = NaN*ones(numTimePoints,numSamples);
asym0p1 = NaN*ones(numTimePoints,numSamples);

%generate numSample trajectories
for iSample = 1 : numSamples

    %generate 100 time point random walk with diffusion coefficient = 1 micron^2/s
    [traj1(:,:,iSample),errFlag] = brownianMotion(2,1,10,0.1);

    %generate 100 time point random walk with diffusion coefficient = 0.1 micron^2/s
    [traj0p1(:,:,iSample),errFlag] = brownianMotion(2,0.1,10,0.1);

end


%calculate the asymmetry parameter for the generated trajectories
%go from length 3 time points to total length in the 1 micron^2/s
%trajectory
for iSample = 1 : numSamples
    for iTP = 3 : numTimePoints
        asym1(iTP,iSample) = asymDetermination(traj1(1:iTP,:,iSample));
    end
end

%calculate the asymmetry parameter for the generated trajectories
%go from length 3 time points to total length in the 0.1 micron^2/s
%trajectory
for iSample = 1 : numSamples
    for iTP = 3 : numTimePoints
        asym0p1(iTP,iSample) = asymDetermination(traj0p1(1:iTP,:,iSample));
    end
end

%get the asymmetry parameter threshold, such that 99% of the sample has a
%lower assymmetry parameter
thresh1 = prctile(asym1,99,2);
thresh0p1 = prctile(asym0p1,99,2);
