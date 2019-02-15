function wrapperMovieExample(movieData)

%% Registration
%Get the indices of any previous threshold processes from this function                                                                              
iProc = movieData.getProcessIndex('ExampleProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    movieData.addProcess(ExampleProcess(movieData));
    iProc = movieData.getProcessIndex('ExampleProcess',1,0);
end

process = movieData.getProcess(iProc);

%% Input check

%% Input/output
p = process.getParameters();
outputFile = fullfile(p.OutputDirectory, 'output.mat');
if ~isdir(p.OutputDirectory), mkdir(p.OutputDirectory), end

% logging input
nChan = numel(movieData.channels_);
inFilePaths = cell(nChan, 1);
inFilePaths{1} = movieData.channels_(1).channelPath_;
process.setInFilePaths(inFilePaths);

% logging output
nChan = numel(movieData.channels_);
outFilePaths = cell(nChan, 1);
outFilePaths{1} = outputFile;
process.setOutFilePaths(outFilePaths);


%% Algorithm

output = zeros(movieData.nFrames_, 1);
for t = 1 : movieData.nFrames_
    I = movieData.getChannel(1).loadImage(t);
    
    % Algorithm
    output(t) = myFunction(I);
        
end
save(outputFile, 'output');

end