app = processC3D_exported;
app.path = 'C:\Users\willi\Desktop\forHannes';
% app.emgLabelCSVPath = 'C:\Users\Willi\SynologyDrive\Hanging_Paper\emg_labels_d2_v2.csv';
fileList = dir(fullfile(app.path, '**', '*.c3d'));

for i = 1: numel(fileList)
    if ~strcmp(fileList(i).name, 'c3dfile.c3d') && isfolder(fullfile(fileList(i).folder, strrep(fileList(i).name, '.c3d', '')))

        pause(0.5);
        app.path = fileList(i).folder;
        copyfile(fullfile(app.path, fileList(i).name), fullfile(app.path, 'tmp.c3d'));

        pause(0.5);
        app.c3dFileName = 'tmp.c3d';
        disp(['processing' fileList(i).name]);

        app.ClimbingCheckBox.Value = 0;
        app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 1;
        app.detectforceplatesautomaticallyCheckBox.Value = 1;
        app.filterorderSpinner.Value = 4;
        app.cutofffrequencySpinner.Value = 20;

        app.filtermarkersCheckBox.Value = 0;
        app.filterGRFsCheckBox.Value = 1;
        app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 1;
        app.filterEMGandexporttostoCheckBox.Value = 0;

        app.processC3Dfile();
        pause(0.5);

        app.CreatefilesButtonPushed(0);

        
        if isfile(fullfile(app.path, 'tmp', 'grf.mot'))
            copyfile(fullfile(app.path, 'tmp', 'grf.mot'), fullfile(app.path, strrep(fileList(i).name, '.c3d', ''), 'grf.mot'));
            rmdir(fullfile(app.path, 'tmp'), 's');
        end
        delete(fullfile(app.path, 'tmp.c3d'));
        s = dir(fullfile(app.path, strrep(fileList(i).name, '.c3d', ''), 'grf.mot'));
        if s.bytes < 100000
            disp('might be wrong');
        end
    end
end
