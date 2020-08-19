function fileList = getAllFiles(dirName)
% get files recursively from subfolders
% taken from:
% http://www.codemiles.com/matlab-examples/get-all-files-in-folder-t7354.html
% for linux turn slash in line 7 !!!
  dirData = dir(dirName);%# Get the data for the current directory
  dirWav = dir(fullfile(dirName, '/*.WAV'));
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirWav.name}';  %'# Get a list of the files
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles(nextDir)];  %# Recursively call getAllFiles
  end
end