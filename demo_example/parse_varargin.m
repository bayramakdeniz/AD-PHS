function varargout = parse_varargin(vargs,defaults)

% parse_varargin.m
% 
% Purpose: parse a varargin cell array and apply defaults where values are
%           not specified by varargin
%
% Usage: process_varargin(vargs,defaults)
%
% Created: 02-22-11 by Tyler Seibert
% See Change Log at end of file for modifications
% 
% Required Arguments:
%   vargs           =   [cell 1xn] varargin cell array from function calling
%                       process_varargin.m
%   defaults    =   [cell mx2] cell array with variable/value pairs
% 

% check defaults and get varnames:
if size(defaults,2) ~= 2
    error('%s: ''defaults'' should be an nx2 cell array with rows corresponding to variable/default pairs\n',mfilename);
else
    varnames = {defaults{:,1}}; % variable names
    defaultvals = {defaults{:,2}}; % default values
end

% check for obvious errors:
if numel(vargs) == 0
    % nothing. This means no optional arguments were provided
elseif rem(numel(vargs),2)
    error('%s: optional arguments should come in pairs (key,value)\n',mfilename);
elseif numel(vargs)/2 > numel(varnames)
    fprintf('%s: more optional arguments than possibilities in ''varnames''\n',mfilename);
    error('%s: %d pairs in varargin, but only %d varnames given\n',mfilename,numel(vargs)/2,numel(varnames));
elseif numel(varnames) ~= nargout
    fprintf('%s: number of outputs must match number of variable names in ''varnames''\n',mfilename);
    error('%s: %d elements in nargout, but %d varnames\n',mfilename,nargout,numel(varnames));
end

% assign all key/value pairs as variables:
for i=1:2:numel(vargs)
    if ~ismember(vargs{i},varnames)
        error('%s: varargin key ''%s'' not found among varnames of ''defaults''\n',mfilename,vargs{i});
    else
        eval([vargs{i} '= vargs{i+1};']);
    end
end

% now assign defaults for all leftovers:
for k=1:numel(varnames)
    name = varnames{k};
    val = defaultvals{k};
    if ~exist(name,'var')
        %fprintf('%s: %s was not assigned. Using default value\n',mfilename,name);
        eval([name '= val;'])
    end
end

% now assign outputs:
for k=1:nargout
    varargout(k) = {eval(varnames{k})};
end

% CHANGE LOG:
% Created 02-22-11 by Tyler Seibert
% 04-05-11 --> changed variable 'v' to 'vargs' to avoid conflict with
%               variables called 'v' in other functions
% 
