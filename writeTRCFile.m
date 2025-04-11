function writeTRCFile(filename, markerData, time)
% writeTRCFile - Export marker data to a .trc file for OpenSim
%
% filename   - string, name of the .trc file to create
% markerData - struct with fields as marker names, each containing an Nx3 matrix
%              representing XYZ positions over time
% time       - Nx1 vector of time stamps (in seconds)

% Get marker names and number of frames
markerNames = fieldnames(markerData);
nMarkers = length(markerNames);
nFrames = length(time);
frameRate = 1 / (time(2) - time(1));

% Check all marker data dimensions
for i = 1:nMarkers
    marker = markerNames{i};
    if size(markerData.(marker), 1) ~= nFrames || size(markerData.(marker), 2) ~= 3
        error(['Marker ' marker ' must be an Nx3 matrix.']);
    end
end

% TRC Header
fileID = fopen(filename, 'w');
fprintf(fileID, 'PathFileType\t4\t(X/Y/Z)\t%s\n', filename);
fprintf(fileID, 'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fileID, '%d\t%d\t%d\t%d\tmm\t%d\t1\t%d\n', frameRate, frameRate, nFrames, nMarkers, frameRate, nFrames);

% Header labels
fprintf(fileID, 'Frame#\tTime');
for i = 1:nMarkers
    fprintf(fileID, '\t%s\t\t', markerNames{i});
end
fprintf(fileID, '\n');

% Sub-header (X/Y/Z)
fprintf(fileID, '\t');
for i = 1:nMarkers
    fprintf(fileID, '\tX%d\tY%d\tZ%d', i, i, i);
end
fprintf(fileID, '\n');

% Write frame data
for frame = 1:nFrames
    fprintf(fileID, '%d\t%.5f', frame, time(frame));
    for i = 1:nMarkers
        coords = markerData.(markerNames{i})(frame, :);
        fprintf(fileID, '\t%.3f\t%.3f\t%.3f', coords(1), coords(2), coords(3));
    end
    fprintf(fileID, '\n');
end

fclose(fileID);
fprintf('TRC file written to %s\n', filename);
end
