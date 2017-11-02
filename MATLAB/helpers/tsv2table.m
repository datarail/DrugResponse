function t = tsv2table(filename,varargin)
% TSV2TABLE Create a table by reading from a tsv file.
%   Use the TSV2TABLE function to create a table by reading column-oriented data
%   from a file.  TSV2TABLE assume a text file with tab delimitation.
%
%   T = TSV2TABLE(FILENAME) creates a table by reading from the tsv file FILENAME.
%   T = TSV2TABLE(FILENAME, 'KeepQuote') creates a table by reading from the tsv
%       file FILENAME. Entries that have quotes (") are kept (by default: deleted)
%   T = TSV2TABLE(FILENAME,..., varargin) creates a table by reading from the tsv
%       file FILENAME. varargin are default READTSV options (see below).
%
%          Reading from a delimited text file creates one variable in T for
%          each column in the file.  Variable names are taken from the first row
%          of the file.  By default, the variables created are either double,
%          if the entire column is numeric, or cell array of strings, if any
%          element in a column is not numeric.  READTABLE converts empty fields
%          in the file to either NaN (for a numeric variable) or the empty string
%          (for a string-valued variable). Insignificant whitespace in the file
%          is ignored.
%
%          Use the following optional parameter name/value pairs to control how
%          data are read from a delimited text file:
%
%
%          'ReadVariableNames'  A logical value that specifies whether or not the
%                          first row (after skipping HeaderRows) of the file is
%                          treated as variable names.  Default is true.
%
%          'ReadRowNames'  A logical value that specifies whether or not the
%                          first column of the file is treated as row names.
%                          Default is false.  If the 'ReadVariableNames' and
%                          'ReadRowNames' parameter values are both true, the
%                          name in the first column of the first row is saved
%                          as the first dimension name for the table.
%
%          'TreatAsEmpty'  One or more strings to be treated as the empty string
%                          in a numeric column.  This may be a character string,
%                          or a cell array of strings.  Table elements
%                          corresponding to these are set to NaN.  'TreatAsEmpty'
%                          only applies to numeric columns in the file, and
%                          numeric literals such as '-99' are not accepted.
%
%          'HeaderLines'   The number of lines to skip at the beginning of the
%                          file.
%
%          'Format'        A format string to define the columns in the file, as
%                          accepted by the TEXTSCAN function.  If you specify 'Format',
%                          you may also specify any of the parameter name/value pairs
%                          accepted by the TEXTSCAN function.  Type "help textscan" for
%                          information about format strings and additional parameters.
%                          Specifying the format can significantly improve speed for
%                          some large files.

if ~isempty(varargin) && (any(strcmpi(varargin,'KeepQuote') | strcmpi(varargin,'KpQ')))
    Rmq = false;
    varargin = varargin(~strcmp(varargin,'KeepQuote') & ~strcmpi(varargin,'KpQ'));
else
    Rmq = true;
end

varargin = [{'Delimiter', '\t', 'FileType', 'text'} varargin];

warning('off','MATLAB:codetools:ModifiedVarnames')
warning('off', 'MATLAB:table:ModifiedVarnames')
t = table.readFromFile(filename,varargin);
warning('on','MATLAB:codetools:ModifiedVarnames')
warning('on', 'MATLAB:table:ModifiedVarnames')

if Rmq
    for i=1:size(t,2)
        if iscellstr(t.(i))
            idx = cellfun(@length, t.(i))>=2;
            t.(i)(idx) = cellfun(@(x)x((x~='"') | [false true(1,length(x)-2) false]),...
                t.(i)(idx),'uniformoutput',0);
        end
    end
end
