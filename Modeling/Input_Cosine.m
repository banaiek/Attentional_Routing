function cosineMatrix = Input_Cosine(Ntotal,samplingRate,duration,frequency,coherence,freqShift)
% Define parameters
frequency = frequency+freqShift*rand(Ntotal,1); % frequency of the cosine wave in Hz

% Calculate the number of samples
numSamples = duration * samplingRate+1;

% Generate the time vector
t = linspace(0, duration, numSamples);

% Initialize the matrix
cosineMatrix = zeros(Ntotal, numSamples);

% Define the range for phase shift based on coherence
% For perfect coherence (coherence = 1), the range is 0
% For no coherence (coherence = 0), the range is 2*pi (full circle)
phaseRange = 2 * pi * (1 - coherence);

% Populate each row of the matrix with a cosine wave with a random phase

    phaseShift = phaseRange * rand(Ntotal,1) ; % random phase shift within the range
    cosineMatrix = cos(2 * pi * frequency .* t + phaseShift);


% The matrix 'cosineMatrix' now contains 1000 cosine waves with controlled phase coherence
imagesc(cosineMatrix)
end