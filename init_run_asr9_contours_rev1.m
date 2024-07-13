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

rev_folder='C:\Local Matlab Data\2.7GHz-FAA\Rev105'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now running the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parallel_flag=0
tf_server_status=0
tf_recalculate=0
%%%%%%%%wrapper_bugsplat_rev1(app,rev_folder,parallel_flag)
%%%%wrapper_bugsplat_rev6(app,rev_folder,parallel_flag,tf_server_status,tf_recalculate)
wrapper_bugsplat_rev7_kml_output(app,rev_folder,parallel_flag,tf_server_status,tf_recalculate)





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







