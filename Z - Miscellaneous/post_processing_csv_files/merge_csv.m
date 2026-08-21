%% Load the files 
clear all
close all

% Load file A and B
[fileA,pathA] = uigetfile({'*.csv'},'Select the localization file A','F:\Louise\DATA\DATA2023');
[fileB,pathB] = uigetfile({'*.csv'},'Select the localization file B',pathA);
filenameA=fullfile(pathA,fileA);
TA = readtable(filenameA,'PreserveVariableNames',true);
A = table2array(TA);
filenameB=fullfile(pathB,fileB);
TB = readtable(filenameB,'PreserveVariableNames',true);
B = table2array(TB);


%% Modify some columns 

B_modified = B;
B_modified(:,1) = B(:,1)+A(end,1); % change frame number 

%% Create new table
Tmerge = [A; B_modified]; 
VarNames = TB.Properties.VariableNames;
%% Add some columns 
% Tmerge = [Tmerge ones(size(Tmerge,1),1)];
% Tmerge(1:size(A,1),end) = zeros(size(A,1),1);

% VarNames{end+1} = 'id_table'

%% Save 
[filename, pathname]=uiputfile({'*.txt' },'Save merged data',pathA);
file_name_full=fullfile(pathname,[filename(1:end-4),'.txt']);
fid=fopen(file_name_full,'w');
display('Saving data .txt')
[r,s]=size(Tmerge);
fprintf(fid,['%g' repmat(['\t %g'],1,s-1) '\n'],Tmerge');
fclose(fid);
display('done')
file_name_full=fullfile(pathname,[filename(1:end-4),'.csv']);
display('Saving data .csv')
Tmerge_table = array2table(Tmerge,'VariableNames',VarNames);
writetable(Tmerge_table,file_name_full);
display('done')