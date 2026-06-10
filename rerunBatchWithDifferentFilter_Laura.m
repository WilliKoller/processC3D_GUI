app = processC3D_exported;
app.path = 'C:\insidebone\00_DATA\Motion_Data\TD31\TD31_2025-11-12_treadmill';
app.emgLabelCSVPath = 'C:\insidebone\00_DATA\Motion_Data\emg_labels.csv';
fileList = dir(fullfile(app.path, '**', '*.c3d'));

processTrials = GetSubDirsFirstLevelOnly(app.path);

for i = 1: numel(fileList)
    if contains(fileList(i).name, 'slow04') && ~strcmp(fileList(i).name, 'c3dfile.c3d') && ~strcmp(fileList(i).name, 'tmp.c3d') && any(contains(processTrials, strrep(fileList(i).name, '.c3d', '')))
%
        % pause(0.5);
        app.path = fileList(i).folder;
        copyfile(fullfile(app.path, fileList(i).name), fullfile(app.path, 'tmp.c3d'));

        % pause(0.5);
        app.c3dFileName = 'tmp.c3d';
        disp(['processing' fileList(i).name]);

        app.AMTITandemTreadmillCheckBox.Value = 1;
        app.filterEMGCheckBox.Value = 1;
        app.renameEMGlabelsCheckBox.Value = 1;

        app.processC3Dfile();
        % pause(0.2);

        app.CreatefilesButtonPushed(0);

        % pause(0.2);

        tmpSteps = GetSubDirsFirstLevelOnly(fullfile(app.path, 'tmp_*'));

        for j = 1 : numel(tmpSteps)
        
            if isfile(fullfile(app.path, tmpSteps{j}, 'EMG_filtered.sto'))
                if isfolder(fullfile(app.path, [strrep(fileList(i).name, '.c3d', '')  tmpSteps{j}(4:end)]))
                    copyfile(fullfile(app.path, tmpSteps{j}, 'EMG_filtered.sto'), fullfile(app.path, [strrep(fileList(i).name, '.c3d', '')  tmpSteps{j}(4:end)], 'EMG_filtered.sto'));
                end
            end

            % % delete(fullfile(app.path, 'tmp.c3d'));
            % s = dir(fullfile(app.path, strrep(fileList(i).name, '.c3d', ''), 'grf.mot'));
            % if s.bytes < 100000
            %     disp('might be wrong');
            % end
            % 
            % if isfile(fullfile(app.path, 'tmp', 'marker_experimental.trc'))
            %     copyfile(fullfile(app.path, 'tmp', 'marker_experimental.trc'), fullfile(app.path, strrep(fileList(i).name, '.c3d', ''), 'marker_experimental.trc'));
            %     % rmdir(fullfile(app.path, 'tmp'), 's');
            % end
            % delete(fullfile(app.path, 'tmp.c3d'));
            % s = dir(fullfile(app.path, strrep(fileList(i).name, '.c3d', ''), 'marker_experimental.trc'));
            % if s.bytes < 100000
            %     disp('might be wrong');
            % end
            
            rmdir(fullfile(app.path, tmpSteps{j}));
        end
        
        delete(fullfile(app.path, 'tmp.c3d'));
    end
end
