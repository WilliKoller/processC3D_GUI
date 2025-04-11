folder = 'C:\Users\willi\ucloud\PhD\GaitModifications\TestData\static01';
filename = 'marker_experimental.trc';
markerData = load_marker_trc(fullfile(folder, 'marker_experimental.trc'));
markerNames = fieldnames(markerData);

% markersOfInterest = {'RASI', 'LASI', 'SACR', 'LKNE', 'RKNE', 'LKNM', 'RKNM', 'LANK', 'RANK', 'LANM', 'RANM'};
markersOfInterest = {'RASI', 'LASI', 'LPSI', 'RPSI', 'LKNE', 'RKNE', 'LKNM', 'RKNM', 'LANK', 'RANK', 'LMED', 'RMED'};
markerDataOfInterest = [];
for i = 1 :numel(markersOfInterest)
    markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterest{i} '_X']));
    markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterest{i} '_Y']));
    markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterest{i} '_Z']));
end

LKJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LKNE'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LKNM'), :, :)))/2;
RKJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RKNE'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RKNM'), :, :)))/2;
LAJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LANK'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LMED'), :, :)))/2;
RAJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RANK'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RMED'), :, :)))/2;
JointCenters.LKJC=LKJC;
JointCenters.RKJC=RKJC;

% add AJC_offset marker based on Bruening et al., 2008
% right AJC_offset
SL=RKJC-RAJC;
SL_distance=sqrt( SL(1,1)*SL(1,1)+SL(1,2)*SL(1,2)+SL(1,3)*SL(1,3) );
offset = 0.027*SL_distance;
for t=1:size(RKJC, 1)
    originRTibia=RAJC;
    [e1TibiaR,e2TibiaR,e3TibiaR]=segmentorientation_r(RKJC(t,:)-RAJC(t,:),squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RANK'), t, :))'-RAJC(t,:));
    rotationmatrix_1=[e1TibiaR;e2TibiaR;e3TibiaR];
    RAJC_offset_rTibia(t,:)=[-offset 0 0];
    RAJC_offset(t,:)=(mldivide(rotationmatrix_1,(RAJC_offset_rTibia(t,:)'))+originRTibia(t,:)')';
    JointCenters.RAJC(t,:)=RAJC_offset(t,:);
end

% left AJC_offset
SL_left=LKJC-LAJC;
SL_left_distance=sqrt( SL_left(1,1)*SL_left(1,1)+SL_left(1,2)*SL_left(1,2)+SL_left(1,3)*SL_left(1,3) );
offset_left = 0.027*SL_left_distance;
for t=1:size(LKJC, 1)
    originLTibia=LAJC;
    [e1TibiaL,e2TibiaL,e3TibiaL]=segmentorientation_l(LKJC(t,:)-LAJC(t,:),squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LANK'), t, :))'-LAJC(t,:));
    rotationmatrix_2=[e1TibiaL;e2TibiaL;e3TibiaL];
    LAJC_offset_lTibia(t,:)=[-offset_left 0 0];
    LAJC_offset(t,:)=(mldivide(rotationmatrix_2,(LAJC_offset_lTibia(t,:)'))+originLTibia(t,:)')';
    JointCenters.LAJC(t,:)=LAJC_offset(t,:);
end

% add Harrington HJC
ASISvector=squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LASI'), :, :))-squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RASI'), :, :));
ASISdistance=sqrt((ASISvector(1,1))^2+(ASISvector(1,2))^2+(ASISvector(1,3))^2);
MidASIS = (squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RASI'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LASI'), :, :)))/2;
if any(strcmp(markersOfInterest, 'SACR'))
    MidPSIS = squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'SACR'), :, :));
else
    MidPSIS = (squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RPSI'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LPSI'), :, :)))/2;
end

PW=ASISdistance;
%   Define the pelvis origin
originPelvis=MidASIS;
for t=1:size(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RASI'), :, :)), 1)
    % HJC definition based on modified Harrington equations (global
    % coordinates) Sangeux, 2015
    H_ap=-0.138*PW-10.4;
    H_v=-0.305*PW-10.9;
    H_ml=0.33*PW+7.3;
    H_ml_left=-H_ml;
    RHJC_Harringtion(t,:) = horzcat(H_ap,H_v,H_ml);
    LHJC_Harringtion(t,:) = horzcat(H_ap,H_v,H_ml_left);

    %   Calculate unit vectors of Pelvis segment relative to global
    [e1Pelvis,e2Pelvis,e3Pelvis]=segmentorientation_1Frame(originPelvis(t,:)-MidPSIS(t,:),(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RASI'), t, :))-squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LASI'), t, :)))');
    rotationmatrix=[e1Pelvis;e2Pelvis;e3Pelvis];
    RHJC(t,:)=(mldivide(rotationmatrix,(RHJC_Harringtion(t,:)'))+originPelvis(t,:)')';
    LHJC(t,:)=(mldivide(rotationmatrix,(LHJC_Harringtion(t,:)'))+originPelvis(t,:)')';
    JointCenters.RHJC(t,:)=RHJC(t,:);
    JointCenters.LHJC(t,:)=LHJC(t,:);
    rotationmatrix=[];
end


uniqueMarkerNames = [];
for m = 1 : numel(markerNames)
    if ~strcmpi(markerNames{m}, 'time') && ~strcmpi(markerNames{m}, 'frame')
        uniqueMarkerNames{end+1} = markerNames{m}(1:end-2);
    end
end
uniqueMarkerNames = unique(uniqueMarkerNames);

markerDataInMat = [];
for m = 1 : numel(uniqueMarkerNames)
    markerDataInMat.(uniqueMarkerNames{m}) = [cell2mat(markerData.([uniqueMarkerNames{m} '_X'])), cell2mat(markerData.([uniqueMarkerNames{m} '_Y'])), cell2mat(markerData.([uniqueMarkerNames{m} '_Z']))];
end

jointCenterNames = fieldnames(JointCenters);
for j = 1 : numel(jointCenterNames)
    markerDataInMat.(jointCenterNames{j}) = JointCenters.(jointCenterNames{j});
end

if isfield(markerData, 'Time')
    time = cell2mat(markerData.Time);
else
    time = cell2mat(markerData.time);
end

writeTRCFile(fullfile(folder, strrep(filename, '.trc', '_withJointCenters.trc')), markerDataInMat, time);
