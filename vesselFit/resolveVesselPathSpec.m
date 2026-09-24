function [fldRoot, node, idx] = resolveVesselPathSpec(vesselV, pathSpec)
% RESOLVEVESSELPATHSPEC  Generic dot-path resolver into ONE scalar vessel struct, consolidating what
% used to be two separate showVessel.m arguments (fld, roi) into a single unified addressing scheme --
% see showVessel.m's own file header PATH SPEC note for the calling convention this serves, and
% drawGaussianPeakContour.m for the companion drawing primitive a resolved 'fit' node feeds into.
%
% SYNTAX: '<fld>.<name1>[(N)].<name2>[(N)]. ...' -- each dot-separated segment is either a plain
% struct field name (optionally followed by a literal '(N)' to index into that field, when it holds a
% struct array) -- see parseDotPathSeg.m for the exact per-segment parsing (promoted from here once
% plotVessels.m's own resolveVesselMetric needed the identical '<name>(N)' form for a numeric row
% index into a .vec matrix) -- OR, when the CURRENT node is already a struct array with a '.label' field and the
% segment text doesn't match any of its own field names, a LABEL MATCH against that array (node.label
% == segment text) instead of a field-name step -- this is what lets 'tsImCentered.rois.excludeManualSinus'
% and 'tsImCentered.rois(1)' both resolve (to the same entry, if that roi happens to be rois(1)):
% the former via label match, the latter via explicit index. The full path down to the actual
% struct-array field also still works ('...patchVesselsSimultaneous.fits(2)') -- but the LAST segment
% may also stop one level short, naming a "wrapper" struct instead (e.g.
% '...patchVesselsSimultaneous(2)', or bare '...patchVesselsSimultaneous') -- see AUTO-DESCEND below.
%
% AUTO-DESCEND (2026-07-24, LAST segment only). If the LAST segment resolves to a SCALAR struct that
% itself has EXACTLY ONE struct-valued field (e.g. fitPatchVessels.m's own outFld wrapper -- .fits is
% its only struct-valued field, alongside plain char/numeric ones like .fitMode/.resnorm), that ONE
% field is silently stepped into BEFORE this segment's own trailing index/defaulting is applied -- so
% '...patchVesselsSimultaneous(2)' (or bare '...patchVesselsSimultaneous') resolves EXACTLY like
% '...patchVesselsSimultaneous.fits(2)' (or '...patchVesselsSimultaneous.fits'), no need to spell out
% '.fits' explicitly. PURELY STRUCTURAL -- no hardcoded field name ('fits' or otherwise): a wrapper
% with zero or 2+ struct-valued fields is left alone (ambiguous/not applicable), requiring the caller to
% spell out the full path in that case -- this resolver still needs NO built-in knowledge of which
% container shapes exist in this codebase (rois vs. fitPatchVessels' own .fits, or any future sibling),
% it just also tolerates stopping one segment short of the array itself when that's unambiguous. Applies
% ONLY to the LAST segment -- an EXPLICIT intermediate '.fits' segment (e.g. the full form above) is
% resolved exactly as before, never auto-descended a second time.
%
% DEFAULTING: if the walk ends (after AUTO-DESCEND, if it applied) on a struct-array field with NO
% trailing index/label (e.g. '...patchVesselsSimultaneous.fits' with nothing after it, or bare
% '...rois'), node defaults to element 1 -- the MAIN/first entry (fitPatchVessels.m's own MAIN peak
% convention for .fits; "first roi" for .rois) -- so a caller wanting the common case never needs to
% spell out '(1)' explicitly. Applies REGARDLESS of how many elements that field happens to hold for
% THIS particular vessel (fixed 2026-07-24 -- previously gated on numel(node)>1, which left idx=[] for a
% struct-array field that happened to hold exactly one element, e.g. a single-peak vessel's own bare
% '.fits' reference; a caller like showVessel.m's own FIT-type idx==1 gold/red gate then treated that as
% "not the main peak" since idx==1 is false when idx=[], coloring a genuinely-main, genuinely-only peak
% red instead of gold -- INCONSISTENTLY across vessels, since a multi-peak vessel's SAME bare reference
% already defaulted correctly. Never applies to the bare-fldRoot case (no segments at all past the root)
% -- see OUTPUT idx below. UNCHANGED by MULTI-INDEX below -- a bare trailing segment (no parens at all)
% still means "element 1", never "every element"; a caller wanting every element spells out '(:)'
% explicitly (plotVessels.m's own vesselPatchFrame route rewrites a bare '...footprints' overlay spec
% to '...footprints(:)' BEFORE calling this function, as ITS OWN opt-in default -- see that file's own
% drawVesselPatchFrameOverlays -- rather than this generic resolver changing what "bare" means for
% every caller/field, which would silently change existing '.rois'/'.fits' bare-reference behavior
% everywhere else this function is already used).
%
% MULTI-INDEX (2026-07-24, Seb's own explicit ask -- LAST segment only, wrapping N>1 struct-array
% entries into ONE overlay call instead of one plotVessels call per index): the LAST segment's own
% '(...)' suffix additionally accepts ':' (every element), an 'end'-relative range ('1:end'), a
% bracketed literal ('[1 3 5]'), or a bare comma list ('1,2,3') -- not just the single bare integer
% parseDotPathSeg.m already supports. See PARSEDOTPATHSEGLAST/RESOLVEIDXEXPR below for the exact
% grammar/evaluation. node/idx (see OUTPUT below) become a struct ARRAY / numeric row vector,
% element-for-element, instead of a scalar/single-int, whenever more than one index is requested --
% every EXISTING single-int/bare call site is unaffected (node/idx stay scalar exactly as before),
% since this is purely additive new '(...)' syntax on the LAST segment. An INTERMEDIATE segment (not
% the last) still only accepts the ORIGINAL single bare-integer form -- multi-selecting partway
% through a walk (continuing to step INTO more than one struct at once) isn't a case this resolver's
% simple field-name-stepping walk supports, so that syntax is reserved for the terminal container only.
%
% INPUT
%   vesselV  : ONE scalar vessel struct (not an array/cell -- callers resolve per-vessel).
%   pathSpec : char, see SYNTAX above. The first segment (before any '.') is always a plain
%              vessel.(fld) field name -- see OUTPUT fldRoot below.
%
% OUTPUT
%   fldRoot : char -- the path's own first segment, i.e. which vessel.(fld) this spec is rooted at
%             (the background image-track field, in showVessel.m's own usage).
%   node    : whatever the path resolved to -- either vessel.(fldRoot) itself (a bare fldRoot with no
%             further segments), a SINGLE struct element (after index/label/default-to-1 narrowing,
%             per DEFAULTING above), or a struct ARRAY (per MULTI-INDEX above, when the last segment's
%             own index expression selected more than one element).
%   idx     : the numeric index(es) actually used for the FINAL struct-array narrowing step (explicit,
%             label-matched, defaulted-to-1, or multi-selected) -- [] if node was never narrowed from a
%             struct array (the bare-fldRoot case); a scalar for every pre-MULTI-INDEX case (unchanged);
%             a row vector, element-for-element matching node, when MULTI-INDEX selected >1 element.
%             Lets a caller distinguish "peak 1/main" from "peak 2" etc. without re-deriving it.
%
%   [fldRoot, node]      = resolveVesselPathSpec(vessel(3), 'tsImCentered.rois(1)');
%   [fldRoot, node]      = resolveVesselPathSpec(vessel(3), 'tsImCentered.rois.excludeManualSinus');
%   [fldRoot, node, idx] = resolveVesselPathSpec(vessel(3), 'tsImCentered.patchVesselsSimultaneous.fits');      % idx=1, main peak
%   [fldRoot, node, idx] = resolveVesselPathSpec(vessel(3), 'tsImCentered.patchVesselsSimultaneous.fits(2)');   % idx=2, 1st secondary
%   [fldRoot, node, idx] = resolveVesselPathSpec(vessel(3), 'tsImCentered.patchVesselsSimultaneous');           % AUTO-DESCEND -- same as .fits above, idx=1
%   [fldRoot, node, idx] = resolveVesselPathSpec(vessel(3), 'tsImCentered.patchVesselsSimultaneous(2)');        % AUTO-DESCEND -- same as .fits(2) above, idx=2
%   [fldRoot, node]      = resolveVesselPathSpec(vessel(3), 'tsImCentered');                                    % bare -- node = vessel(3).tsImCentered itself
%   [fldRoot, node, idx] = resolveVesselPathSpec(vessel(3), 'tsIm.footprints(:)');                              % MULTI-INDEX -- node/idx cover EVERY footprint
%   [fldRoot, node, idx] = resolveVesselPathSpec(vessel(3), 'tsIm.footprints(1,3)');                            % MULTI-INDEX -- node/idx = footprints([1 3])

    assert(isstruct(vesselV) && isscalar(vesselV), 'resolveVesselPathSpec:notScalarVessel', ...
        'vesselV must be a single scalar vessel struct -- resolve per-vessel, not per-array/cell.');
    segs = strsplit(pathSpec, '.');
    [name0, idx0] = parseDotPathSeg(segs{1});
    assert(isempty(idx0), 'resolveVesselPathSpec:indexedRoot', ...
        'the root segment (''%s'') names a vessel.(fld) field and cannot carry an index.', segs{1});
    fldRoot = name0;
    assert(isfield(vesselV, fldRoot), 'resolveVesselPathSpec:noRootField', ...
        'vessel has no field ''%s''.', fldRoot);
    node = vesselV.(fldRoot);
    idx = [];

    for i = 2:numel(segs)
        if i==numel(segs)
            [name, idxExpr] = parseDotPathSegLast(segs{i});   % see file header MULTI-INDEX
        else
            [name, idxExpr] = parseDotPathSeg(segs{i});   % UNCHANGED -- single bare-integer only
        end
        if isstruct(node) && isfield(node, name)
            node = node.(name);
            idx = [];
            if i==numel(segs)
                node = autoDescendWrapper(node);   % see file header AUTO-DESCEND -- LAST segment only
            end
            if ~isempty(idxExpr)
                segIdx = resolveIdxExpr(idxExpr, numel(node));
                assert(all(segIdx>=1) && all(segIdx<=numel(node)), 'resolveVesselPathSpec:indexOutOfRange', ...
                    '''%s(%s)'' out of range (1..%d) in ''%s''.', name, mat2str(segIdx), numel(node), pathSpec);
                node = node(segIdx);
                idx = segIdx;
            end
        elseif isstruct(node) && isfield(node,'label')
            assert(isempty(idxExpr), 'resolveVesselPathSpec:labelWithIndex', ...
                'segment ''%s'' looks like a label match (not a field of the current node) but also carries an index, in ''%s''.', segs{i}, pathSpec);
            k = find(strcmp({node.label}, name), 1);
            assert(~isempty(k), 'resolveVesselPathSpec:noLabelMatch', ...
                'no entry with .label==''%s'' found (%d entries) resolving ''%s''.', name, numel(node), pathSpec);
            node = node(k);
            idx = k;
        else
            error('resolveVesselPathSpec:badSegment', ...
                'segment ''%s'' does not resolve against the current node (class %s) in ''%s''.', segs{i}, class(node), pathSpec);
        end
    end

    % DEFAULTING (see file header): ended on an un-narrowed struct (array or scalar) -- take element 1.
    % Guarded on numel(segs)>1 -- the bare-fldRoot case (no segments past the root at all) must keep
    % idx=[] (see OUTPUT idx above), not default to 1.
    if numel(segs)>1 && isstruct(node) && isempty(idx)
        node = node(1);
        idx = 1;
    end
end

% ---------------------------------------------------------------------------
function node = autoDescendWrapper(node)
    % See file header AUTO-DESCEND. A no-op unless node is a SCALAR struct with EXACTLY ONE
    % struct-valued field, in which case node steps into that one field.
    if ~isstruct(node) || numel(node)~=1; return; end
    fn = fieldnames(node);
    if isempty(fn); return; end
    isStructField = structfun(@isstruct, node);
    structFields = fn(isStructField);
    if isscalar(structFields)
        node = node.(structFields{1});
    end
end

% ---------------------------------------------------------------------------
function [name, idxExpr] = parseDotPathSegLast(seg)
    % Same '<name>(N)' syntax as parseDotPathSeg.m, for the LAST path segment ONLY -- additionally
    % accepts a MULTI-element index expression inside the parens (see file header MULTI-INDEX): ':'
    % (every element), an 'end'-relative range ('1:end'), a bracketed literal ('[1 3 5]'), or a bare
    % comma list ('1,2,3'). Returns the RAW, unevaluated inner text (resolved into concrete indices by
    % RESOLVEIDXEXPR below, once the target array's own length is known) rather than parsing it here --
    % this function only splits <name> from its own '(...)' suffix, capturing the WHOLE parenthesized
    % suffix as ONE group and manually stripping the parens afterward (tok{2}(2:end-1)) -- the SAME
    % technique parseDotPathSeg.m itself uses, and for the SAME reason: empirically confirmed (not
    % merely assumed -- see that file's own header) that a capturing group NESTED inside another
    % capturing/non-capturing group is silently dropped by MATLAB's own regexp 'tokens'+'once', so the
    % inner index text can't be captured directly as its own nested group here either.
    %
    % INPUT
    %   seg : one (LAST) dot-path segment, e.g. 'footprints(:)', 'rois(1,2)', 'fits([1 3])', 'tsIm'.
    %
    % OUTPUT
    %   name    : seg's own leading token, without any '(...)' suffix.
    %   idxExpr : '' if seg carried no '(...)' suffix at all (bare -- caller's own DEFAULTING applies,
    %             UNCHANGED from parseDotPathSeg.m's own bare-segment behavior) -- otherwise the raw
    %             text INSIDE the parens, unevaluated (see RESOLVEIDXEXPR).
    %
    %   [name, idxExpr] = parseDotPathSegLast('footprints(:)');   % name='footprints', idxExpr=':'
    %   [name, idxExpr] = parseDotPathSegLast('rois(1,2,3)');     % name='rois', idxExpr='1,2,3'
    %   [name, idxExpr] = parseDotPathSegLast('fits(2)');         % name='fits', idxExpr='2'
    %   [name, idxExpr] = parseDotPathSegLast('tsIm');            % name='tsIm', idxExpr=''
    %
    % See also parseDotPathSeg, resolveIdxExpr.
    tok = regexp(seg, '^(\w+)(\(.*\))?$', 'tokens', 'once');
    assert(~isempty(tok), 'resolveVesselPathSpec:badSegmentSyntax', 'malformed path segment ''%s''.', seg);
    name = tok{1};
    if numel(tok)<2 || isempty(tok{2})
        idxExpr = '';
    else
        idxExpr = tok{2}(2:end-1);   % strip the captured '(' and ')'
    end
end

% ---------------------------------------------------------------------------
function segIdx = resolveIdxExpr(idxExpr, n)
    % Evaluates a PARSEDOTPATHSEGLAST-returned index-expression into a concrete positive-integer row
    % vector, against n (the target struct array's own numel -- what a bare ':' or an 'end'-relative
    % range resolve against). A NUMERIC idxExpr (parseDotPathSeg.m's own single-int-or-[] return, for
    % every INTERMEDIATE segment) passes straight through unchanged -- only a CHAR idxExpr (the LAST
    % segment's own extended grammar) is actually evaluated here. str2num (NOT eval directly, though
    % str2num itself uses eval internally per its own docs) evaluates the remaining MATLAB-numeric-
    % literal syntax (a range, a bracketed list, or a bare comma list, auto-wrapped in brackets) once
    % ':' and 'end' are substituted -- acceptable here since pathSpec strings are caller-authored MATLAB
    % source in this codebase already (e.g. plotVessels.m's own 'fname(args...)' colormap grammar
    % already accepts arbitrary caller-supplied argument text the same way), not externally-supplied
    % untrusted input.
    %
    % INPUT
    %   idxExpr : numeric (pass-through) OR raw index-expression text (PARSEDOTPATHSEGLAST's own
    %             idxExpr, non-empty) -- ':', a comma list ('1,2,3'), a bracketed list ('[1 3]'), or a
    %             range possibly using the literal word 'end' ('1:end').
    %   n       : the target struct array's own numel -- what ':'/'end' resolve against.
    %
    % OUTPUT
    %   segIdx : positive-integer row vector (1-based).
    %
    %   segIdx = resolveIdxExpr(':', 5);        % [1 2 3 4 5]
    %   segIdx = resolveIdxExpr('1:end', 5);    % [1 2 3 4 5]
    %   segIdx = resolveIdxExpr('1,2,3', 5);    % [1 2 3]
    %   segIdx = resolveIdxExpr('[1 3]', 5);    % [1 3]
    %   segIdx = resolveIdxExpr(3, 5);          % 3 (numeric pass-through)
    if isnumeric(idxExpr)
        segIdx = idxExpr;
        return
    end
    inner = strtrim(idxExpr);
    if strcmp(inner, ':')
        segIdx = 1:n;
        return
    end
    inner = regexprep(inner, '\<end\>', num2str(n));
    if isempty(regexp(inner, '^\[.*\]$', 'once'))
        inner = ['[' inner ']'];
    end
    segIdx = str2num(inner); %#ok<ST2NM>
    assert(isnumeric(segIdx) && isvector(segIdx) && ~isempty(segIdx) && all(segIdx==round(segIdx)), ...
        'resolveVesselPathSpec:badIndexExpr', 'could not parse index expression ''(%s)''.', idxExpr);
    segIdx = segIdx(:)';
end
