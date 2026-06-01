function writeTRCFile(filename, markerData, time, xyz)
% writeTRCFile - Export marker data to a .trc file for OpenSim
%
% filename   - string, name of the .trc file to create
% markerData - struct with fields as marker names, each containing an Nx3 matrix
%              representing XYZ positions over time
% time       - Nx1 vector of time stamps (in seconds)
% xyz        - defines the structure of markerData. 

if ~isempty(xyz) && xyz == 1
    markerNames = fieldnames(markerData);
    uniqueMarkerNames = [];
    for m = 1 : numel(markerNames)
        if ~strcmpi(markerNames{m}, 'time') && ~strcmpi(markerNames{m}, 'frame')
            uniqueMarkerNames{end+1} = markerNames{m}(1:end-2);
        end
    end
    uniqueMarkerNames = unique(uniqueMarkerNames);

    markerDataInMat = [];
    for m = 1 : numel(uniqueMarkerNames)
        try
            markerDataInMat.(uniqueMarkerNames{m}) = [cell2mat(markerData.([uniqueMarkerNames{m} '_X'])), cell2mat(markerData.([uniqueMarkerNames{m} '_Y'])), cell2mat(markerData.([uniqueMarkerNames{m} '_Z']))];
            if ischar(markerDataInMat.(uniqueMarkerNames{m}))
                rows = size(markerDataInMat.(uniqueMarkerNames{m}), 1);
                markerDataInMat.(uniqueMarkerNames{m}) = [];
                markerDataInMat.(uniqueMarkerNames{m})(1:rows, 1:3) = NaN;
            end
        catch
            try
                markerDataInMat.(uniqueMarkerNames{m}) = [markerData.([uniqueMarkerNames{m} '_X']), markerData.([uniqueMarkerNames{m} '_Y']), markerData.([uniqueMarkerNames{m} '_Z'])];
            end
        end
    end

    markerData = markerDataInMat;
    markerNames = fieldnames(markerData);
end

% Get marker names and number of frames
markerNames = fieldnames(markerData);
nMarkers = length(markerNames);
nFrames = length(time);
frameRate = 1 / (time(2) - time(1));

% Check all marker data dimensions
maxTravellingDistance = 0;
for i = 1:nMarkers
    marker = markerNames{i};
    if size(markerData.(marker), 1) ~= nFrames || size(markerData.(marker), 2) ~= 3
        error(['Marker ' marker ' must be an Nx3 matrix.']);
    else
        if isnumeric(markerData.(marker))
            maxTravellingDistance = max(maxTravellingDistance, pdist2(markerData.(marker)(1, :), markerData.(marker)(end, :)));
        end
    end
end

% TRC Header
fileID = fopen(filename, 'w');
fprintf(fileID, 'PathFileType\t4\t(X/Y/Z)\t%s\n', filename);
fprintf(fileID, 'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
if maxTravellingDistance < 20
    fprintf(fileID, '%d\t%d\t%d\t%d\tm\t%d\t1\t%d\n', frameRate, frameRate, nFrames, nMarkers, frameRate, nFrames);
else
    fprintf(fileID, '%d\t%d\t%d\t%d\tmm\t%d\t1\t%d\n', frameRate, frameRate, nFrames, nMarkers, frameRate, nFrames);
end

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
        if isnumeric(coords)
            fprintf(fileID, '\t%.3f\t%.3f\t%.3f', coords(1), coords(2), coords(3));
        else
            coords = cell2mat(coords);
            if isnumeric(coords)
                fprintf(fileID, '\t%.3f\t%.3f\t%.3f', coords(1), coords(2), coords(3));
            else
                fprintf(fileID, '\t%.3f\t%.3f\t%.3f', nan, nan, nan);
            end
        end
    end
    fprintf(fileID, '\n');
end

fclose(fileID);
fprintf('TRC file written to %s\n', filename);
end
