function fh = openfh(path, varargin)
%OPENFH Open file.
%
%   FH = OPENFH(PATH)
%   FH = OPENFH(PATH, PERMISSIONS)
%   FH = OPENFH(PATH, PERMISSIONS, MACHINEFORMAT)
%   FH = OPENFH(PATH, PERMISSIONS, MACHINEFORMAT, ENCODING)
%
%   OPENFH if a light wrapper for the built-in FOPEN function.  See
%   the documentation for FOPEN for a description of its arguments.
%
%   The main difference between OPENFH and FOPEN is that, while FOPEN
%   returns -1 on error, OPENFH fails outright.
%
%   Also, whereas FOPEN can be called with a file descriptor as its
%   first and only argument, the first argument to OPENFH must be a
%   string.

    narginchk(1, 3);
    if ~isstr_(path); error('first argument must be a string'); end
    [fh, msg] = fopen(path, varargin{:});
    if fh == -1; error(msg); end
end
