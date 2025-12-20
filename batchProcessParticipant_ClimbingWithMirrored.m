app = processC3D_exported;
participantIDs = [14 15 24 26 ];

for p = 1 : numel(participantIDs)
    if participantIDs(p) < 10
        app.path = ['C:\Users\Willi\SynologyDrive\Climbing_Campus\Data\C0' num2str(participantIDs(p))];
    else
        app.path = ['C:\Users\Willi\SynologyDrive\Climbing_Campus\Data\C' num2str(participantIDs(p))];
    end

    % for first partcipant, EMG labels are different!!!
    % app.path = ['C:\Users\Willi\SynologyDrive\Climbing_Campus\Data\C' num2str(participantID)];
    app.emgLabelCSVPath = 'C:\Users\Willi\SynologyDrive\Hanging_Paper\emg_labels_d2_v2.csv';
    fileList = dir([app.path '\*.c3d']);
    %%
    for i = 1 : numel(fileList)
        app.c3dFileName = fileList(i).name;
        disp(['processing ' app.c3dFileName]);

        app.ClimbingCheckBox.Value = 1;
        app.invertforcesensorsCheckBox.Value = 1;
        app.mirrorCheckBox.Value = 0;
        app.rotatearoundyaxisDropDown.Value = '180°';

        app.detectwalkingdirectionautomaticallybymarkerCheckBox.Value = 0;
        app.detectforceplatesautomaticallyCheckBox.Value = 0;
        app.filtermarkersCheckBox.Value = 0;
        app.filterGRFsCheckBox.Value = 0;
        app.setGRFstozerooutsideregionoffootcontactCheckBox.Value = 0;

        app.filterEMGandexporttostoCheckBox.Value = 1;
        app.normalizeEMGCheckBox.Value = 0;

        % CHANGE
        app.renameEMGlabelsCheckBox.Value = 0;

        try
            app.processC3Dfile();

            app.CreatefilesButtonPushed(0);

            app.processC3Dfile();

            app.mirrorCheckBox.Value = 1;

            app.CreatefilesButtonPushed(0);
        catch e
            disp(['an error occured in ' app.c3dFileName ' ' e.message])
        end
    end
end


% forceIDs = zeros(size(forces, 1), 1);
% for i = 1 :size(forces, 1)
%    if forces(i, 2) > 0
%        % if norm(forces(i, :)) < 30
%        %     forceIDs(i) = 1;
%        % end
%    else
%        forceIDs(i) = 0;
%    end
% end