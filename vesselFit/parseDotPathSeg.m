function [name, idx] = parseDotPathSeg(seg)
% PARSEDOTPATHSEG  Parses ONE dot-path segment into its plain name and optional trailing '(N)' index
% -- promoted from resolveVesselPathSpec.m's own local parseSeg once a SECOND caller (plotVessels.m's
% resolveVesselMetric, for a literal numeric '.vec(N)' row index) needed the identical '<name>(N)'
% parsing -- see this codebase's established promote-once-2+-callers-need-it convention
% (getNestedField.m, resolveIrfCore.m, etc.).
%
% SYNTAX: '<name>' or '<name>(N)' -- see resolveVesselPathSpec.m's own file header SYNTAX for the
% calling convention this parses one segment of.
%
% INPUT
%   seg : one dot-path segment, e.g. 'rois(2)', 'vec', 'tsRSF'.
%
% OUTPUT
%   name : seg's own leading token, without any '(N)' suffix.
%   idx  : the parsed integer N, or [] if seg carried no '(N)' suffix.
%
%   [name, idx] = parseDotPathSeg('rois(2)');   % name='rois', idx=2
%   [name, idx] = parseDotPathSeg('vec');       % name='vec',  idx=[]
%
% See also resolveVesselPathSpec, plotVessels.

    % The index suffix is its OWN capturing group, INCLUDING its parens -- '(\(\d+\))?', not
    % '(?:\((\d+)\))?' -- confirmed empirically (not merely assumed) that MATLAB's regexp 'tokens'
    % does NOT surface a capturing group nested inside a non-capturing '(?:...)' one at all (the token
    % silently vanishes rather than coming back empty), so nesting the digit-only group inside a
    % non-capturing paren wrapper as originally written here silently discarded every index -- caught
    % by a real test (resolveVesselPathSpec('...rois(2)') was returning rois(1), not rois(2)), not by
    % code review alone. Parens are stripped below instead. Also: 'tokens'+'once' omits a
    % non-participating group's OWN slot entirely (numel(tok)==1, not 2-with-an-empty-placeholder) when
    % the suffix is absent, so a missing 2nd element means "no index", same as an empty one would.
    tok = regexp(seg, '^(\w+)(\(\d+\))?$', 'tokens', 'once');
    assert(~isempty(tok), 'parseDotPathSeg:badSegmentSyntax', 'malformed path segment ''%s''.', seg);
    name = tok{1};
    if numel(tok)<2 || isempty(tok{2})
        idx = [];
    else
        idx = str2double(tok{2}(2:end-1));   % strip the captured '(' and ')'
    end
end
