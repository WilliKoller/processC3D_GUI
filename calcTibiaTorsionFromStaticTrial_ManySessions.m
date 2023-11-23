clear;
% baseFolder = 'C:\Users\willi\ucloud\DATA_Hans_PhD';
% additionalSubDir = 0;
baseFolder = 'C:\Users\willi\Documents\UniDataLokal\DEMAND\VICON_DATA_SCHMELZ';
additionalSubDir = 1;

participantFolders = GetSubDirsFirstLevelOnly(baseFolder);
for p = 1 : numel(participantFolders)
    sessions = [];
    if additionalSubDir == 1
        sessions = GetSubDirsFirstLevelOnly(fullfile(baseFolder, participantFolders{p}));
    else
        sessions{1} = '';
    end
    for s = 1 : numel(sessions)
        if isempty(sessions{s})
            sessionFolder = fullfile(baseFolder, participantFolders{p});
        else
            sessionFolder = fullfile(baseFolder, participantFolders{p}, sessions{s});
        end
        trials = lower(GetSubDirsFirstLevelOnly(sessionFolder));
        staticTrial = trials(or(or(or(contains(trials, 'cal'), contains(trials, 'stand')), contains(trials, 'static')), contains(trials, 'initial')));
        
        if ~isempty(staticTrial)
            disp(['Processing ' participantFolders{p} ' - ' staticTrial{1}]);
            folder = fullfile(sessionFolder, staticTrial{1});

            markerData = load_marker_trc(fullfile(folder, 'marker_experimental.trc'));
            markerNames = fieldnames(markerData);

            markersOfInterest = {'RASI', 'LASI', 'LKNE', 'RKNE', 'LKNM', 'RKNM', 'LANK', 'RANK', 'LANM', 'RANM', 'SACR'};
            markersOfInterestAlternative = {'RASI', 'LASI', 'LKNE', 'RKNE', 'LMKNE', 'RMKNE', 'LANK', 'RANK', 'LMMA', 'RMMA', 'SACR'};
            markersOfInterestAlternative2 = {'RASI', 'LASI', 'LKNE', 'RKNE', 'LKNEM', 'RKNEM', 'LANK', 'RANK', 'LANKM', 'RANKM', 'SACR'};
            markerDataOfInterest = [];
            for i = 1 :numel(markersOfInterest)
                if strcmp(markersOfInterest{i}, 'SACR')
                    try
                        markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterest{i} '_X']));
                        markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterest{i} '_Y']));
                        markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterest{i} '_Z']));
                    catch
                        markersOfInterest{i} = 'RPSI';
                        markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterest{i} '_X']));
                        markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterest{i} '_Y']));
                        markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterest{i} '_Z']));
                        i = i + 1;
                        markersOfInterest{i} = 'LPSI';
                        markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterest{i} '_X']));
                        markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterest{i} '_Y']));
                        markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterest{i} '_Z']));
                    end
                else
                    try
                        markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterest{i} '_X']));
                        markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterest{i} '_Y']));
                        markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterest{i} '_Z']));
                    catch
                        try
                            markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterestAlternative{i} '_X']));
                            markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterestAlternative{i} '_Y']));
                            markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterestAlternative{i} '_Z']));
                        catch
                            markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterestAlternative2{i} '_X']));
                            markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterestAlternative2{i} '_Y']));
                            markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterestAlternative2{i} '_Z']));
                        end
                    end
                end
            end

            LKJC=(squeeze(markerDataOfInterest(contains(markersOfInterest, 'LKNE'), :, :))+squeeze(markerDataOfInterest(contains(markersOfInterest, 'LKNM'), :, :)))/2;
            RKJC=(squeeze(markerDataOfInterest(contains(markersOfInterest, 'RKNE'), :, :))+squeeze(markerDataOfInterest(contains(markersOfInterest, 'RKNM'), :, :)))/2;
            LAJC=(squeeze(markerDataOfInterest(contains(markersOfInterest, 'LANK'), :, :))+squeeze(markerDataOfInterest(contains(markersOfInterest, 'LANM'), :, :)))/2;
            RAJC=(squeeze(markerDataOfInterest(contains(markersOfInterest, 'RANK'), :, :))+squeeze(markerDataOfInterest(contains(markersOfInterest, 'RANM'), :, :)))/2;
            JointCenters.LKJC=LKJC;
            JointCenters.RKJC=RKJC;

            %% calc tibia torsion 
            GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0; ...
                norm(cross(A,B)) dot(A,B)  0; ...
                0              0           1];

            FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

            UU = @(Fi,G) Fi*G*inv(Fi);

            LkneeAxis = mean(squeeze(markerDataOfInterest(contains(markersOfInterest, 'LKNE'), :, :)) - squeeze(markerDataOfInterest(contains(markersOfInterest, 'LKNM'), :, :)), 1);
            LkneeAxis = LkneeAxis / norm(LkneeAxis);
            RkneeAxis = mean(squeeze(markerDataOfInterest(contains(markersOfInterest, 'RKNE'), :, :)) - squeeze(markerDataOfInterest(contains(markersOfInterest, 'RKNM'), :, :)), 1);
            RkneeAxis = RkneeAxis / norm(RkneeAxis);

            LankleAxis = mean(squeeze(markerDataOfInterest(contains(markersOfInterest, 'LANK'), :, :)) - squeeze(markerDataOfInterest(contains(markersOfInterest, 'LANM'), :, :)), 1);
            LankleAxis = LankleAxis / norm(LankleAxis);
            RankleAxis = mean(squeeze(markerDataOfInterest(contains(markersOfInterest, 'RANK'), :, :)) - squeeze(markerDataOfInterest(contains(markersOfInterest, 'RANM'), :, :)), 1);
            RankleAxis = RankleAxis / norm(RankleAxis);

            RtibiaLongitundialAxis = mean(RKJC - RAJC, 1);
            LtibiaLongitundialAxis = mean(LKJC - LAJC, 1);

            % rotate axis in tibia coordinate system
            % right tibia
            a=RtibiaLongitundialAxis';
            a = a /norm(a);
            b=[0 1 0]';
            U = UU(FFi(a,b), GG(a,b));

            RkneeAxisRotated = U * RkneeAxis';
            RankleAxisRotated = U * RankleAxis';

            u = RkneeAxisRotated([1 3]);
            v = RankleAxisRotated([1 3]);
            RtibiaTorsion = atan2d(u(1)*v(2)-u(2)*v(1),u(1)*v(1)+u(2)*v(2));

            % left tibia
            a=LtibiaLongitundialAxis';
            a = a /norm(a);
            b=[0 1 0]';
            U = UU(FFi(a,b), GG(a,b));

            LkneeAxisRotated = U * LkneeAxis';
            LankleAxisRotated = U * LankleAxis';

            u = LkneeAxisRotated([1 3]);
            v = LankleAxisRotated([1 3]);
            LtibiaTorsion = atan2d(u(1)*v(2)-u(2)*v(1),u(1)*v(1)+u(2)*v(2)) * -1;

            % positve value = external rotation
            % negative value = internal rotation
            save(fullfile(folder, 'tibiaTorsion.mat'), 'RtibiaTorsion', "LtibiaTorsion");
        end
    end
end