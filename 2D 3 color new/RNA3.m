%%%%%%%% Latest version 3 color RNA quantification %%%%%%%%%%
% Start with a folder and get a list of all subfolders.
% Finds and prints names of all PNG, JPG, and TIF images in 
% that folder and all of its subfolders.
clc;% Clear the command window
clear; 
workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;
pixelsize = 39.682539;
radius = 300; %radius of exclusion (for the same channel)
dist = 1; %radius of inclusion (for two different channels)
alchk = 0; % 1 if you want to check alignment, 0 if you don't
mask1 = 'Cymask.tif'; %cytoplasmidc mask file name
mask2 = 'Nucmask.tif'; % use '' in case there is only one mask
mrnafile = ["5p.loc"; "mid.loc"; "3p.loc"];
reference=["mid"];

% Define a starting folder.
% start_path = fullfile(matlabroot, '\toolbox\images\imdemos');
% Ask user to confirm or change.
% topLevelFolder = uigetdir(start_path);
% if topLevelFolder == 0
% 	return;
% end
% Get list of all subfolders.
topLevelFolder = pwd; %current working directory
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

if numberOfFolders == 1
    cprintf('err','error: No subfolders\n');
    return
end
jst =[];
% baseFileNames = [];
% Process all image files in those folders.
for k = 1 : numberOfFolders-1
	% Get this folder and print it out.
	thisFolder = listOfFolderNames{k+1};
% 	fprintf('Processing folder %s\n', thisFolder);
	cd(thisFolder)
%     thisFolders = ls;
%    
%     for lm = 1:size(thisFolders,1)
%         jst{lm} = char2cell(thisFolders(lm,:));
%     end
  if exist('3p.loc', 'file') == 2 && exist('Cymask.tif','file') == 2 
   RNA_coloc3(mask1, mask2, mrnafile, pixelsize, radius, dist, alchk,reference);
  else
      cprintf('err','NO files to process\n')
  end
        cd ..
% 	% Get PNG files.
% 	filePattern = sprintf('%s/*.png', thisFolder);
% 	baseFileNames = dir(filePattern);
% 	% Add on TIF files.
% 	filePattern = sprintf('%s/*.tif', thisFolder);
% 	baseFileNames = [baseFileNames; dir(filePattern)];
% 	% Add on JPG files.
% 	filePattern = sprintf('%s/*.jpg', thisFolder);
% 	baseFileNames = [baseFileNames; dir(filePattern)];
% 	numberOfImageFiles = length(baseFileNames);
% 	% Now we have a list of all files in this folder.
% 	
% 	if numberOfImageFiles >= 1
% 		% Go through all those image files.
% 		for f = 1 : numberOfImageFiles
% 			fullFileName = fullfile(thisFolder, baseFileNames(f).name);
% 			fprintf('     Processing image file %s\n', fullFileName);
% 		end
% 	else
% 		fprintf('     Folder %s has no image files in it.\n', thisFolder);
% 	end
end

cprintf('Comments','Done\n')
