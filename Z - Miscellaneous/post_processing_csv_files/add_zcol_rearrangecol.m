% Load file 
[file,path] = uigetfile({'*.csv'},'Select the localization file','F:\Louise\DATA\DATA2023');
filename=fullfile(path,file);
T = readtable(filename,'PreserveVariableNames',true);
A = table2array(T);

%% 4 col: x-y-z-t
T_new = ones(size(A,1),4);
T_new(:,1:2) = A(:,2:3);
T_new(:,4) = A(:,1);
T_new_titles =  {'x [nm]' 'y [nm]' 'z [nm]' 'frame'};
%% Save 
[filename, pathname]=uiputfile({'*.txt' },'Save new table',path);
file_name_full=fullfile(pathname,[filename(1:end-4),'.txt']);
fid=fopen(file_name_full,'w');
display('Saving data .txt')
[r,s]=size(T_new);
fprintf(fid,['%g' repmat(['\t %g'],1,s-1) '\n'],T_new');
fclose(fid)
display('done')
file_name_full=fullfile(pathname,[filename(1:end-4),'.csv']);
display('Saving data .csv')
T_table = array2table(T_new,'VariableNames',T_new_titles);
writetable(T_table,file_name_full);
display('done')