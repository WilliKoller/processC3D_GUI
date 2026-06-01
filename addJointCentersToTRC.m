folder = 'I:\data_speiseng\P01\pre\Static03_rot';

markerData = load_marker_trc(fullfile(folder, 'marker_experimental.trc'));
markerNames = fieldnames(markerData);

% markersOfInterest = {'RASI', 'LASI', 'SACR', 'LKNE', 'RKNE', 'LKNM', 'RKNM', 'LANK', 'RANK', 'LANM', 'RANM'};
markersOfInterest = {'RASI', 'LASI', 'LPSI', 'RPSI', 'LKNE', 'RKNE', 'LKNEM', 'RKNEM', 'LANK', 'RANK', 'LANKM', 'RANKM'};
markerDataOfInterest = [];
for i = 1 :numel(markersOfInterest)
    markerDataOfInterest(i, :, 1) = cell2mat(markerData.([markersOfInterest{i} '_X']));
    markerDataOfInterest(i, :, 2) = cell2mat(markerData.([markersOfInterest{i} '_Y']));
    markerDataOfInterest(i, :, 3) = cell2mat(markerData.([markersOfInterest{i} '_Z']));
end

LKJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LKNE'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LKNM'), :, :)))/2;
RKJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RKNE'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RKNM'), :, :)))/2;
LAJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LANK'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'LANM'), :, :)))/2;
RAJC=(squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RANK'), :, :))+squeeze(markerDataOfInterest(strcmp(markersOfInterest, 'RANM'), :, :)))/2;
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

% remove /1000 if m or mm mismatch!
markerData.LHJC_X = JointCenters.LHJC(:, 1) / 1000;
markerData.LHJC_Y = JointCenters.LHJC(:, 2) / 1000;
markerData.LHJC_Z = JointCenters.LHJC(:, 3) / 1000;
markerData.RHJC_X = JointCenters.RHJC(:, 1) / 1000;
markerData.RHJC_Y = JointCenters.RHJC(:, 2) / 1000;
markerData.RHJC_Z = JointCenters.RHJC(:, 3) / 1000;

markerData.LKJC_X = JointCenters.LKJC(:, 1) / 1000;
markerData.LKJC_Y = JointCenters.LKJC(:, 2) / 1000;
markerData.LKJC_Z = JointCenters.LKJC(:, 3) / 1000;
markerData.RKJC_X = JointCenters.RKJC(:, 1) / 1000;
markerData.RKJC_Y = JointCenters.RKJC(:, 2) / 1000;
markerData.RKJC_Z = JointCenters.RKJC(:, 3) / 1000;

markerData.LAJC_X = JointCenters.LAJC(:, 1) / 1000;
markerData.LAJC_Y = JointCenters.LAJC(:, 2) / 1000;
markerData.LAJC_Z = JointCenters.LAJC(:, 3) / 1000;
markerData.RAJC_X = JointCenters.RAJC(:, 1) / 1000;
markerData.RAJC_Y = JointCenters.RAJC(:, 2) / 1000;
markerData.RAJC_Z = JointCenters.RAJC(:, 3) / 1000;

markerDataTime = cell2mat(markerData.Time);
writeTRCFile(fullfile(folder, 'marker_experimental_with_JointCenters.trc'), markerData, markerDataTime, 1);