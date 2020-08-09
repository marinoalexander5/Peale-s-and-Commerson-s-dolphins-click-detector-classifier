function [dnums] = wavname2date( filename )
% Adapted from Triton [http://cetus.ucsd.edu/technologies_Software.html]
%
% Parses .wav file names and converts them to matlab datenums
% Works on single or multiple file names
% Supports 4 filename formats:
%   1.  yymmdd-HHMMSS
%   2.  yymmdd_HHMMSS
%   3.  yyyymmdd_HHMMSS 
%   4.  yymmddHHMMSS
% ouptuts:
%   dnums: for time manipulation - getting detection time relative to date
%    
% check https://strftime.org/ for formating in panads when importing
% for now transforming all cases to avisoft format for easier handling
% date_str = regexprep(date_str, '[\_.-]', '');
% pd.to_datetime(str(date_str), format='%y%m%d%H%M%S') pandas conversion

% regexp(fname,'\d{4}[-_]\d{4}','match','split')

% start with the default delimiter

date_str = regexp(filename,'\d{6}[-]\d{6}','match');

if isempty(date_str) % using underscores presumably
    date_fmt = 'yymmdd_HHMMSS';
    date_str = regexp(filename,'\d{6}[_]\d{6}','match');
end

if isempty(date_str) % not just and underscore problem, try PAMGuard filename
    date_fmt = 'yyyymmdd_HHMMSS'; % PAMGuard default file format 
    date_str = regexp(filename,'\d{8}[_]\d{6}','match');
%     if ~isempty(date_str)
%         fprintf('Using PAMGuard filename format yyyymmdd_HHMMSS');
%     end
end

if isempty(date_str) % not a PAMguard file, try avisoft/SoundTrap filename
    date_fmt = 'yymmddHHMMSS';
    date_str = regexp(filename,'\d{12}','match' );
%     if ~isempty(date_str)
%         fprintf('Using avisoft filename format yymmddHHMMSS');
%     end
end

if isempty(date_str)
    fprintf('Unknown filename date format.  Please use one of the following:');
    fprintf('*yymmdd-HHMMSS*.wav\n');
    fprintf('*yymmdd_HHMMSS*.wav\n');
    fprintf('*yyyymmdd_HHMMSS*.wav\n');
    fprintf('*yymmddHHMMSS*.wav\n'); 
    date_fmt = 'yymmdd-HHMMSS';
end 

dnums = datenum(date_str, date_fmt);
%% smarter to do all at once, need to figure out when
% dnums = cellfun(@(x)datenum(x,date_fmt),date_strs,'UniformOutput',false);
% dnums = cell2mat(dnums)'; % output of cellfun is a cell
