clear;
clc;
close all;
close all force;
app=NaN(1);  %%%%%%%%%This is to allow for Matlab Application integration.
format shortG
format longG
top_start_clock=clock;
folder1='C:\Local Matlab Data\2.7GHz-FAA';
cd(folder1)
pause(0.1)
addpath('C:\Local Matlab Data\General_Terrestrial_Pathloss')  %%%%%%%%This is where we put the other github repositories
addpath('C:\Local Matlab Data\Generic_Bugsplat')
addpath('C:\Local Matlab Data\GMF_functions')



'Coordination zones for ASR-8/9/11-Nexrad'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Propagation Inputs
FreqMHz=2700; %%%%%%%%MHz
reliability=10;
array_reliability_check=reliability;
confidence=50;
Tpol=1; %%%polarization for ITM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load('BW100_FullSpectrum_Test1.mat')
%%%2640MHz Center
bs_tf_mask=vertcat(horzcat(2640,0),horzcat(2687,0.1),horzcat(2689,5),horzcat(2690,19),horzcat(2692,41),horzcat(2695,49),horzcat(2697,54),horzcat(2699,61),horzcat(2705,90),horzcat(2770,98));
%bs_measure=horzcat(FreqMHz',CorrMag');
%[interp_fold_min]=fold_and_min_rev1(app,bs_measure);
% % close all;
% % figure;
% % hold on;
% % plot(FreqMHz,CorrMag,'-ob')
% % plot(bs_tf_mask(:,1),-1*bs_tf_mask(:,2),'--or')
% % %plot(interp_fold_min(:,1),-1*interp_fold_min(:,2),'-g')
% % grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GMF inputs:  'Need to pull in all the ASR-8/9/11 locations and lowest frequency'
gmf_MinMHz = 2700;
gmf_MaxMHz = 3000;
rev_num=2 
%[gmf_table]=pull_gmf_excel_rev1(app,gmf_MinMHz,gmf_MaxMHz,rev_num);


%%%%%%%%%%%%%%%%Can't connect to the GMF because of the PIV card.
filename_gmf_table='gmf_table_2700_3500_4052024.mat';
load(filename_gmf_table,'gmf_table')

gmf_header=gmf_table.Properties.VariableNames



'First do the ASR-9'
%%%%%%%%%%%%%%%%%%%%%%%%Filter: Find Equipment with ASR-9 and FAA GMF Serial
cell_gmf=table2cell(gmf_table);
serial_col_idx=find(matches(gmf_header,'SER'))
serial_xeq_idx=find(matches(gmf_header,'XEQ'))

%%%%Need to fill in the empty equipment rows
header_str='XEQ'
[cell_gmf]=fix_blank_gmf_string(app,gmf_header,header_str,cell_gmf);

%%%%%%%%%%%Convert the XAH from a string to a number
header_string='XAH';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);
[cell_gmf]=gmf_convert_str2num(app,cell_gmf,gmf_header,header_string); %%%%%Convert str 2 num



%%%%%%%%%%%Find the overlap
faa_row_idx=find(contains(cell_gmf(:,serial_col_idx),'FAA'));
asr9_row_idx=find(contains(cell_gmf(:,serial_xeq_idx),'ASR-9'));
faa_asr9_idx=intersect(faa_row_idx,asr9_row_idx);
%%%cell_gmf(faa_asr9_idx,[serial_col_idx,serial_xeq_idx])


%%%%%%%%Need the Frequency, Lat/Lon, Antenna Height, of all ASR-9
str_header_keep=cell(6,1);
str_header_keep{1}='SER';
str_header_keep{2}='XEQ';
str_header_keep{3}='FRQMHz';
str_header_keep{4}='XLatDD';
str_header_keep{5}='XLonDD';
str_header_keep{6}='XAH';

%%%%%%%Only keep this data: Contour: Second Column in the Cell
num_keep_cols=length(str_header_keep);
keep_col_idx=NaN(num_keep_cols,1);
for col_idx=1:1:num_keep_cols
    keep_col_idx(col_idx)=find(matches(gmf_header,str_header_keep{col_idx}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%This is the ASR-9 information to use.
cell_gmf_asr9=cell_gmf(faa_asr9_idx,keep_col_idx);
table_gmf_asr9=cell2table(cell_gmf_asr9);

%%%%%%%%%%Find ASR-9 that are close to each other and take the lowest frequnecy
array_asr9_latlon=cell2mat(cell_gmf_asr9(:,[4,5]));

[idx_knn]=knnsearch(array_asr9_latlon,array_asr9_latlon,'k',2); %%%Find Nearest Neighbor
size(idx_knn)
knn_group_latlon=array_asr9_latlon(idx_knn(:,2),:);
size(knn_group_latlon)
knn_dist_bound=deg2km(distance(knn_group_latlon(:,1),knn_group_latlon(:,2),array_asr9_latlon(:,1),array_asr9_latlon(:,2)));%%%%Calculate Distance
max_knn_dist=ceil(max(knn_dist_bound))
max(knn_dist_bound)
off_idx=find(knn_dist_bound>0)

%%%%%%%%There are is one set of lat/lon off by 63m, fixing
array_asr9_latlon(off_idx,:)
idx_knn(off_idx,:)
array_asr9_latlon(off_idx(2),:)=array_asr9_latlon(off_idx(1),:);
cell_gmf_asr9(off_idx(2),4)=cell_gmf_asr9(off_idx(1),4);
cell_gmf_asr9(off_idx(2),5)=cell_gmf_asr9(off_idx(1),5);

array_asr9_latlon=round(array_asr9_latlon,4);

[~,ia_idx,ic_idx]=unique(array_asr9_latlon,'rows');
size(ia_idx)

%%%%%%%%%%%%%%Find the lowest frequency per location
num_uni_eut=length(ia_idx); %%%Number
uni_loc_gmf_asr9=cell(num_uni_eut,6);
for i=1:1:num_uni_eut
    lia_idx=find(ismember(array_asr9_latlon,array_asr9_latlon(ia_idx(i),:), 'rows'));
    uni_loc_gmf_asr9(i,:)=cell_gmf_asr9(lia_idx(1),:);
    uni_loc_gmf_asr9{i,3}=round(min(cell2mat(cell_gmf_asr9(lia_idx,3)))); 
end

%%%%uni_loc_gmf_asr9


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate FDR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf_calc_fdr=0%1
rx_freq_mhz=max(cell2mat(uni_loc_gmf_asr9(:,3)))
tx_freq_mhz=2640; %%%%Center Frequency of the Base Stations

filename_fdr=strcat('asr9_','array_fdr_',num2str(rx_freq_mhz),'.mat');
[var_exist_fdr]=persistent_var_exist_with_corruption(app,filename_fdr);

if tf_calc_fdr==1
    var_exist_fdr=0
end

if var_exist_fdr==2
    tic;
    load(filename_fdr,'array_fdr')
    toc;
else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FDR Curves
    %%%%'import the asr-9 curve'
    filename_if_curves='ASR_9_IF.xlsx'
    tic;
    table_ars9_if=readtable(filename_if_curves);
    toc;

    array_asr9_if=table2array(table_ars9_if);
    rx_if_asr9=vertcat(horzcat(0,0),horzcat(0.2,1),horzcat(0.3,2.6),horzcat(0.4,4.6),horzcat(0.5,6.8),horzcat(1,20),horzcat(2,39),horzcat(3,52),horzcat(5,67),horzcat(8,86));

    %%%%[interp_fold_min]=fold_and_min_rev1(app,array_asr9_if);

    close all;
    figure;
    hold on;
    plot(array_asr9_if(:,1),array_asr9_if(:,2),'-r')
    plot(-1*rx_if_asr9(:,1),-1*rx_if_asr9(:,2),'-ob')
    %%%plot(interp_fold_min(:,1),-1*interp_fold_min(:,2),'-g')
    grid on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FDR Inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5GBND418 (Frank): 2590-2690 MHz --> 2640MHz Center
    %%%%%%Brian's update
    array_tx_rf=fliplr(horzcat(0, 35, 50,51, 55, 60, 75, 85)); %%%%Frequency MHz (Base Station) [Half Bandwidth]
    array_tx_mask=fliplr(horzcat(0,1, 19, 42, 50, 70, 95, 100)); %%%%%%%dB Loss
    tx_extrap_loss=-15;

% % % %%%%%%%%%Rev 101-107: Old Mask
% % %     array_tx_rf=fliplr(horzcat(0, 50, 55, 80, 150)); %%%%Frequency MHz (Base Station) [Half Bandwidth]
% % %     array_tx_mask=fliplr(horzcat(0, 0.1, 50, 70, 100)); %%%%%%%dB Loss
    tx_extrap_loss=-20; %%%%%%%%%TX Extrapolation Slope dB/Decade -60dB (This is generous)

    array_rx_if=fliplr(rx_if_asr9(:,1)'); %%%%Frequency MHz (ASR-9) [Need to check these numbers]
    array_rx_loss=fliplr(rx_if_asr9(:,2)'); %%%%%%%dB Loss (ASR-9)

    array_rx_if=fliplr(horzcat(0,.664,1.3,2,4, 5.5, 10, 20))/2; %%%%Frequency MHz (ASR-9) [Need to check these numbers]
    array_rx_loss=fliplr(horzcat(0,3,10,20, 40, 50, 68.75, 100)); %%%%%%%dB Loss (ASR-9)
    rx_extrap_loss=-20; %%%%%%%%%RX Extrapolation Slope dB/Decade 60dB (This is generous)

    
    fdr_freq_separation=abs(tx_freq_mhz-rx_freq_mhz)
    fdr_calc_mhz=ceil(fdr_freq_separation*1.25)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate FDR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [FDR_dB,ED,VD,OTR,DeltaFreq,~,trans_mask_lte]=FDR_ModelII_app(app,fdr_calc_mhz,array_tx_rf,array_rx_if,array_tx_mask,array_rx_loss,tx_extrap_loss,rx_extrap_loss);
    toc;

    zero_idx=nearestpoint_app(app,0,DeltaFreq);
    array_fdr=horzcat(DeltaFreq(zero_idx:end)',FDR_dB(zero_idx:end));

    fdr_idx=nearestpoint_app(app,fdr_freq_separation,array_fdr(:,1));
    fdr_dB=array_fdr(fdr_idx,:)  %%%%%%Frequency, FDR Loss
    save(filename_fdr,'array_fdr')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FDR Plot
    close all;
    figure;
    hold on;
    plot(array_fdr(:,1),array_fdr(:,2),'-b','LineWidth',2,'DisplayName','FDR Loss')
    %xline(fdr_dB(1),'-g','LineWidth',2,'DisplayName','Delta F [MHz]')
    %yline(fdr_dB(2),'-r','LineWidth',2,'DisplayName','FDR Loss [dB]')
    legend('Location','northwest')
    title({strcat('FDR: ASR-9 and Base Station')})
    grid on;
    xlabel('Frequency Offset [MHz]')
    ylabel('FDR [dB]')
    filename1=strcat('FDR1_NewTx.png');
    saveas(gcf,char(filename1))
end


% % % %%%%%%%%%%%%%%%%%% % % Bob's Notes
% % % % % % ASR-9 3 dB BW 0.564 MHz
% % % % % % ASR-9 NF 3 dB
% % % % % % ASR-9 Antenna Gain 32 dBi.
% % % % % % RCVR Noise Power is -113.5 dBm.
% % % % % % Protection Criteria I/N of -6 dB
% % % % % % Maximum Interference power from 5G base station in radar receiver is then -119.3 dBm.
%%%%%%radar_threshold=-119.3;

%%%%%%%%ASR-9: Other Inputs
rx_ant_heigt_m=10; %%%%%%10 meters
rx_nf=3%4.5;  %%%%%%%ASR-9 NF in dB
rx_ant_gain_mb=32; %%%%%%Main Beam gain of ASR-9 in dBi
in_ratio=-10%-6; %%%%%I/N Ratio -6dB
rx_bw_mhz=0.664;%0.564;
[radar_threshold]=rx_threshold_calc_rev1(app,rx_bw_mhz,rx_nf,in_ratio)
dpa_threshold=floor(rx_ant_gain_mb-radar_threshold)  %%%%%%%%%-149dBm
%%%%%-174+10*log10(rx_bw_mhz*10^6)+rx_nf+in_ratio




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Simulation Input Parameters to change
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % rev=101; %%%%%%Example at a Single Location: locations with 2705
% % grid_spacing=2;  %%%%km
% % sim_radius_km=10; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% % bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% % array_bs_eirp_reductions=bs_eirp;
% % bs_height=30; %%%%%Height in m
% % array_mitigation=0:10:20;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% % %loc_idx1=find(contains(uni_loc_gmf_asr9(:,1),'FAA 891029'));
% % array_freq=cell2mat(uni_loc_gmf_asr9(:,3))
% % loc_idx1=find(array_freq<=2705)
% % cell_locations=uni_loc_gmf_asr9([loc_idx1],:)
% % tf_calc_rx_angle=0;
% sim_folder1=folder1
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rev=102; %%%%%%Example at a Single Location: locations with 2705
% grid_spacing=2;  %%%%km
% sim_radius_km=10; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% array_bs_eirp_reductions=bs_eirp;
% bs_height=30; %%%%%Height in m
% array_mitigation=0:10:20;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% %loc_idx1=find(contains(uni_loc_gmf_asr9(:,1),'FAA 891029'));
% array_freq=cell2mat(uni_loc_gmf_asr9(:,3))
% loc_idx1=find(array_freq<=2705)
% cell_locations=uni_loc_gmf_asr9([loc_idx1(1)],:)
% tf_calc_rx_angle=0;
% sim_folder1=folder1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % rev=103; %%%%%%Example at a Single Location: locations with 2705
% % grid_spacing=5;  %%%%km
% % sim_radius_km=40; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% % bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% % array_bs_eirp_reductions=bs_eirp;
% % bs_height=30; %%%%%Height in m
% % array_mitigation=0:20:20;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% % array_freq=cell2mat(uni_loc_gmf_asr9(:,3))
% % loc_idx1=find(array_freq<=2705)
% % cell_locations=uni_loc_gmf_asr9([loc_idx1(1)],:)
% % tf_calc_rx_angle=0;
% sim_folder1=folder1
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rev=104; %%%%%%Example at a Single Location: locations with 2705
% grid_spacing=5;  %%%%km
% sim_radius_km=60; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% array_bs_eirp_reductions=bs_eirp;
% bs_height=30; %%%%%Height in m
% array_mitigation=0:20:20;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% %loc_idx1=find(contains(uni_loc_gmf_asr9(:,1),'FAA 891029'));
% array_freq=cell2mat(uni_loc_gmf_asr9(:,3))
% loc_idx1=find(array_freq<=2705)
% cell_locations=uni_loc_gmf_asr9([loc_idx1(1:2)],:)
% tf_calc_rx_angle=0;
% sim_folder1=folder1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % rev=105; %%%%%%All locations
% % grid_spacing=5;  %%%%km
% % sim_radius_km=70; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% % bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% % array_bs_eirp_reductions=bs_eirp;
% % bs_height=30; %%%%%Height in m
% % array_mitigation=0:20:20;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% % cell_locations=uni_loc_gmf_asr9;
% % tf_calc_rx_angle=0;
% sim_folder1=folder1
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rev=106; %%%%%%Example at a Single Location: locations with 2705: Runtime
% grid_spacing=5;  %%%%km
% sim_radius_km=60; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% array_bs_eirp_reductions=bs_eirp;
% bs_height=30; %%%%%Height in m
% array_mitigation=0:20:60;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% %loc_idx1=find(contains(uni_loc_gmf_asr9(:,1),'FAA 891029'));
% array_freq=cell2mat(uni_loc_gmf_asr9(:,3))
% loc_idx1=find(array_freq<=2705)
% cell_locations=uni_loc_gmf_asr9([loc_idx1(1)],:)
% tf_calc_rx_angle=0;
% sim_folder1=folder1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % rev=107; %%%%%%All locations
% % grid_spacing=1;  %%%%km
% % sim_radius_km=80; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% % bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% % array_bs_eirp_reductions=bs_eirp;
% % bs_height=30; %%%%%Height in m
% % array_mitigation=0:10:60;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% % cell_locations=uni_loc_gmf_asr9;
% % tf_calc_rx_angle=0;
% % sim_folder1='Z:\MATLAB\2.7GHz-FAA'
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % rev=108; %%%%%%All: Debug FAA820962
% % grid_spacing=1;  %%%%km
% % sim_radius_km=80; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% % bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% % array_bs_eirp_reductions=bs_eirp;
% % bs_height=30; %%%%%Height in m
% % array_mitigation=0:10:60;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% % loc_idx1=find(contains(uni_loc_gmf_asr9(:,1),'FAA 820962'));
% % cell_locations=uni_loc_gmf_asr9([loc_idx1(1)],:)
% % tf_calc_rx_angle=0;
% % sim_folder1='Z:\MATLAB\2.7GHz-FAA'
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % rev=109; %%%%%%All locations: New FDR
% % grid_spacing=5;  %%%%km
% % sim_radius_km=80; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
% % bs_eirp=75; %%%%%EIRP [dBm/100MHz]
% % array_bs_eirp_reductions=bs_eirp;
% % bs_height=30; %%%%%Height in m
% % array_mitigation=0:10:60;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
% % array_freq=cell2mat(uni_loc_gmf_asr9(:,3))
% % loc_idx1=find(array_freq<=2705)
% % cell_locations=uni_loc_gmf_asr9([loc_idx1(1)],:)
% % tf_calc_rx_angle=0;
% % sim_folder1=folder1
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rev=110; %%%%%%All locations:: New FDR: Trying to find the max Distance
grid_spacing=1;  %%%%km
sim_radius_km=80; %%%%%%%%Placeholder distance --> Simplification: This is an automated calculation, but requires additional processing time.
bs_eirp=75; %%%%%EIRP [dBm/100MHz]
array_bs_eirp_reductions=bs_eirp;
bs_height=30; %%%%%Height in m
array_mitigation=0:10:60;  %%%%%%%%% in dB: [0, 10, 20, 30, 40, 50]
cell_locations=uni_loc_gmf_asr9;
tf_calc_rx_angle=0;
sim_folder1='Z:\MATLAB\2.7GHz-FAA'
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Link Budget
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Pr = PT + GT + GR -FDR- Pl – System Loss  (System loss is typical value of 2 dB)
required_loss=bs_eirp+dpa_threshold


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Piece together the cell_sim_data
cell_sim_data(:,1)=cell_locations(:,1);

array_latlon=cell2mat(cell_locations(:,[4,5]));
[num_loc,~]=size(array_latlon)
cell_latlon=cell(num_loc,1);
array_freq=cell2mat(cell_locations(:,3))
array_required_pathloss=NaN(num_loc,1);
for i=1:1:num_loc
    cell_latlon{i}=horzcat(array_latlon(i,:),rx_ant_heigt_m);

    freq_separation=abs(array_freq(i)-tx_freq_mhz);
    fdr_idx=nearestpoint_app(app,freq_separation,array_fdr(:,1));
    fdr_dB=array_fdr(fdr_idx,2);  %%%%%%Frequency, FDR Loss
    array_required_pathloss(i)=required_loss-fdr_dB;
end
cell_sim_data(:,2)=cell_latlon; %%%Polygon
cell_sim_data(:,3)=cell_latlon; %%%%Protection Points
cell_sim_data(:,4)=num2cell(array_required_pathloss);

cell_sim_data'



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Create a Rev Folder
% % cd(folder1);
% % pause(0.1)
% % tempfolder=strcat('Rev',num2str(rev));
% % [status,msg,msgID]=mkdir(tempfolder);
% % rev_folder=fullfile(folder1,tempfolder);
% % cd(rev_folder)
% % pause(0.1)

cd(sim_folder1);
pause(0.1)
tempfolder=strcat('Rev',num2str(rev));
[status,msg,msgID]=mkdir(tempfolder);
rev_folder=fullfile(sim_folder1,tempfolder);
cd(rev_folder)
pause(0.1)



tic;
save('reliability.mat','reliability')
save('array_reliability_check.mat','array_reliability_check')
save('confidence.mat','confidence')
save('FreqMHz.mat','FreqMHz')
save('Tpol.mat','Tpol')
save('sim_radius_km.mat','sim_radius_km')
save('grid_spacing.mat','grid_spacing')
save('array_bs_eirp_reductions.mat','array_bs_eirp_reductions')
save('bs_height.mat','bs_height')
save('cell_sim_data.mat','cell_sim_data')
save('array_mitigation.mat','array_mitigation')
save('tf_calc_rx_angle.mat','tf_calc_rx_angle')
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Saving the simulation files in a folder for the option to run from a server
%%%%%%%%%For Loop the Locations
[num_locations,~]=size(cell_sim_data)
base_id_array=1:1:num_locations; %%%%ALL
table([1:num_locations]',cell_sim_data(:,1))

for base_idx=1:1:num_locations
    strcat(num2str(base_idx/num_locations*100),'%')

    temp_single_cell_sim_data=cell_sim_data(base_idx,:);
    data_label1=temp_single_cell_sim_data{1};
    data_label1 = strrep(data_label1,' ','') 

    %%%%%%%%%Make a Folder each Location/System
    cd(rev_folder);
    pause(0.1)
    tempfolder2=strcat(data_label1);
    [status,msg,msgID]=mkdir(tempfolder2);
    sim_folder=fullfile(rev_folder,tempfolder2);
    cd(sim_folder)
    pause(0.1)

    tic;
    base_polygon=temp_single_cell_sim_data{2};
    save(strcat(data_label1,'_base_polygon.mat'),'base_polygon')

    base_protection_pts=temp_single_cell_sim_data{3};
    save(strcat(data_label1,'_base_protection_pts.mat'),'base_protection_pts')

    required_pathloss=temp_single_cell_sim_data{4};
    save(strcat(data_label1,'_required_pathloss.mat'),'required_pathloss')
    toc;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now running the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Go run the simulation'



end_clock=clock;
total_clock=end_clock-top_start_clock;
total_seconds=total_clock(6)+total_clock(5)*60+total_clock(4)*3600+total_clock(3)*86400;
total_mins=total_seconds/60;
total_hours=total_mins/60;
if total_hours>1
    strcat('Total Hours:',num2str(total_hours))
elseif total_mins>1
    strcat('Total Minutes:',num2str(total_mins))
else
    strcat('Total Seconds:',num2str(total_seconds))
end
cd(folder1)
pause(0.1)







