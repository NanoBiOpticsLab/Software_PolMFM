% code to rename images and correct any numbering issue especially so that
% fishquant can identify the correct number of digits to store the frame
% number
% The frame number here is considered to be placed between two "_" and that
% is is the only number in the name with this configuration 

% if different configuration , you need to readjust it 

% written by Bassam HAJJ
% 27/07/2023


clear all
close all


selpath = uigetdir('C:\','Select the folder location');
fileList =  dir(fullfile(selpath, '*.tif'));
l=size(fileList,1); 
DigitsNumber=numel(num2str(l));
%DigitsNumber=3;
%i=1
display('Running ----- Please don''t do anything----- Please wait-----')
for i=1:l
    filename=fileList(i).name ; %'ma_trix_1_corr_1.tif'
    [filepath, f,ext] = fileparts(filename);
    k = strfind(f,'_');
    
    if i==1 %run one check on two file names to check weather this is the frame number or just a number of a cell or field of view
        filename2=fileList(i+1).name ; %take a second name in the list
        [filepath2, f2,ext2] = fileparts(filename2);
        k2 = strfind(f2,'_');
        for j=1: length(k)-1  % here we assume all names are arranged in a similar way
            str=f(k(j)+1:k(j+1)-1);
            str2=f2(k2(j)+1:k2(j+1)-1);
            if  ~strcmp(str,str2) && ~isnan(str2double(str))% check if different numbering and double checking that this is a number not a string
                j_selected=j; 
            end
        end
        
    else 
%             for j=1: length(k)-1 
%                 str=f(k(j)+1:k(j+1)-1);
%                 if ~isnan(str2double(str))%double checking that this is the number of frames
%                     newf = [f(1:k(j)),sprintf(['%0' num2str(DigitsNumber) 'd'],str2num(str)),f(k(j+1):end)] ;
%                     newname=fullfile(fileList(i).folder,[newf,ext]);
%                     oldname=fullfile(fileList(i).folder,fileList(i).name);
%                     if  ~strcmp(newname,oldname)
%                         movefile(oldname, newname); %safety feature in case some files are already with the correct number of digits
%                     end
%                     break  %assuming it is the correct place for the numer
%                 end
%             end
            j= j_selected;
            str=f(k(j)+1:k(j+1)-1);
            newf = [f(1:k(j)),sprintf(['%0' num2str(DigitsNumber) 'd'],str2num(str)),f(k(j+1):end)] ;
            newname=fullfile(fileList(i).folder,[newf,ext]);
            oldname=fullfile(fileList(i).folder,fileList(i).name);
            if  ~strcmp(newname,oldname)
                movefile(oldname, newname); %safety feature in case some files are already with the correct number of digits
            end
             
            
    end
            
end      
        
display('----------DONE------------')        
        
        
        
%     C=mat2cell(filename([k(1:end-1)'+1,k(2:end)'-1]), 1, length(k)-1);
%     
%     @i (filename(k(i)+1:k(i+1)-1))
%     isnan(str2double(str))
    