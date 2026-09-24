function ax = showVessel(vessel, fld, opts)
% SHOWVESSEL  Primitive display: the vessel patch (time-/run-averaged, physical mm) with overlay(s) --
% ROI outline(s) and/or fitPatchVessels.m peak contour(s) -- drawn on top, if any are requested/exist. A
% small, reusable building block other plotting functions can call for a consistent "show the vessel +
% its overlays" rendering -- the actual patch-panel PIXEL-DRAWING mechanics (imagesc/CLim/colormap/
% toX,toY) are themselves factored out into plotVesselPatch.m, which fitVesselDiag.m now also calls
% directly (see that file's own DIAG STRUCT note) rather than reimplementing the same handful of lines a
% second time -- showVessel.m stays the one that additionally knows about overlay LABELS/legends/the
% GRID LAYOUT below, none of which plotVesselPatch.m needs to know about. Deliberately UNOPINIONATED
% about figure lifecycle otherwise: no dark theme, no docking, no save-to-disk (unlike the
% diagnostic-figure family, plotRoiFitFinal.m/fitVesselDiag.m).
% Plots into opts.ax if given (SCALAR vessel only); else a single vessel gets its own new figure+axes,
% and MULTIPLE vessels are tiled into a GRID instead -- see GRID LAYOUT below.
%
% PATH SPEC (2026-07-24) -- fld's own former sibling argument, roi, is GONE: fld itself now carries
% BOTH what fld used to mean AND what roi used to mean, via TWO calling forms:
%   - fld is char (or empty) -- TODAY'S ORIGINAL, unchanged behavior: a plain background field name (or
%     the default-resolved last 'ts*' field when empty), with overlay auto-defaulting to the FIRST (up
%     to) 3 already-present roi labels on that field -- see DEFAULT OVERLAY below. Backward compatible
%     with every existing call site.
%   - fld is a cellstr -- NEW unified addressing: each entry is a full dot-path spec, resolved via
%     resolveVesselPathSpec.m (see that file for the exact SYNTAX/DEFAULTING/AUTO-DESCEND rules), naming
%     either an ROI ('<fld>.rois(N)' or '<fld>.rois.<label>') or a fitPatchVessels.m PEAK
%     ('<fld>.<outFld>.fits(N)', 'N' defaulting to 1 -- the MAIN peak -- when omitted, e.g.
%     '<fld>.<outFld>.fits' alone). The trailing '.fits' may also be OMITTED (2026-07-24,
%     resolveVesselPathSpec.m's own AUTO-DESCEND) -- '<fld>.<outFld>(N)' or bare '<fld>.<outFld>'
%     resolves identically to spelling out '.fits(N)'/'.fits'. EVERY entry must share the same root
%     field (the panel's own background image, asserted) -- one axes shows one background image, so
%     mixing roots within one spec list is a caller mistake, not a supported "combine two flds" feature.
%     NO auto-defaulting in this form -- an explicit cellstr is a fully-specified overlay request; pass
%     {} for "background only, no overlay" (background then falls back to the DEFAULT FLD rule below,
%     since an empty cellstr carries no root of its own to read it from).
%   Examples (see file header EXAMPLE list below for more):
%     showVessel(vessel, 'tsIm');                                                    % old form, unchanged
%     showVessel(vessel, {'tsIm.rois.lumenGauss', 'tsIm.rois.lumenGauss_tissue'});    % old roi-list equivalent
%     showVessel(vessel, {'tsIm.rois(1)', 'tsIm.patchVesselsSimultaneous.fits'});     % roi + MAIN peak
%     showVessel(vessel, {'tsIm.rois(1)', 'tsIm.patchVesselsSimultaneous'});          % same MAIN peak, '.fits' omitted
%     showVessel(vessel, {'tsIm.patchVesselsSimultaneous.fits(2)'});                  % 1st SECONDARY peak only
%     showVessel(vessel, {'tsIm.irf.irfIm.stats.perVoxel'});                          % activation-map overlay, see STAT-TYPE below
%
% STAT-TYPE OVERLAY (2026-08-18) -- a PATH SPEC entry resolving to a struct with .F/.P fields (fitIRF.m's
% own per-voxel whole-block test layout, e.g. '<fld>.irf.irfIm.stats.perVoxel') is drawn as a per-voxel
% ACTIVATION MAP via drawStatMapOverlay.m: colored by .F magnitude (opts.statOverlay.colormap, a
% SEQUENTIAL colormap since an F-statistic has no sign), opaque only where .P < opts.statOverlay.pThresh,
% transparent (grayscale patch shows through) elsewhere. Requires the vessel to have been fit with
% opts.stats.enable=true and 'perVoxel' in opts.stats.tests (fitIRF.m's own default) -- an
% unfit/unavailable path resolves like any other missing overlay spec (warning + skipped, see PATH SPEC
% above), not an error. Inspired by tools/vasomoTools' now-sunset plotActOld.m/addOverlay.m (grayscale
% underlay + colored, significance-thresholded functional map), reimplemented as a single truecolor
% `image` call (see drawStatMapOverlay.m's own RECIPE note for why no second linked axes is needed here).
% NOT part of the ROI color/style cycling (STYLE below) or the DEFAULT OVERLAY auto-list -- opt-in only,
% via an explicit PATH SPEC entry.
%
% GRID LAYOUT (more than one vessel, opts.ax not given): tiled 3 rows x 5 cols (15) per figure,
% OVERFLOWING to a new figure (same 3x5 tiledlayout) every 15 vessels -- a skipped vessel (missing
% image data) still consumes its own slot (left empty), so vessel v always lands in the SAME slot of
% the SAME figure regardless of any skips elsewhere in the array. In the PATH SPEC form, BOTH the
% background field and the overlay list are resolved ONCE (from vessel(1)) and shared by every vessel in
% the grid -- in the char/empty form, fld still resolves INDEPENDENTLY per vessel (unchanged from
% before this file's own PATH SPEC addition) while only the DEFAULT OVERLAY below is forced shared.
%
% INPUT
%   vessel : vessel struct (scalar or array), or a cell array of vessel structs -- same convention as
%            fitVessel.m. Needs vessel.(fld).im/.vSize; vessel.(fld).rois is OPTIONAL (no overlay drawn
%            if missing/empty -- "with an overlay IF ANY EXIST/ARE REQUESTED").
%   fld    : char (background field, DEFAULT OVERLAY below) | cellstr (explicit PATH SPEC list above) |
%            empty (default fld -- the LAST field on vessel(v) whose name starts with 'ts', resolved
%            independently per vessel, even in GRID LAYOUT -- same default as fitVessel.m -- plus
%            DEFAULT OVERLAY below). A bare char that itself contains a '.' (e.g.
%            'tsIm.rois.lumenManual', typed without the {}) is auto-promoted to a single-entry
%            cellstr PATH SPEC (2026-08-17) -- a literal vessel.(fld) field name can never contain a
%            dot, so a dotted char is unambiguously a PATH SPEC typed without its wrapping braces,
%            not a genuine (if unusual) field name; silently resolving it as ONE background field
%            named e.g. 'tsIm.rois.lumenManual' (which never exists) previously just skipped every
%            vessel with a "no <fld> data" message, no error, easy to miss (Seb's own report,
%            2026-08-17: "I would have expected line 389 to 391 of doIt_human.m to just work").
%   opts   : struct --
%       .ax    target axes to draw into -- ONLY valid when vessel is scalar (errors otherwise, since
%              one axes can't host more than one vessel). Default: a new figure+axes (scalar vessel)
%              or the GRID LAYOUT above (multiple vessels).
%       .color RGB triplet | Nx3 matrix | cell array of RGB triplets -- the color PALETTE cycled across
%              ROI-type overlay entries ONLY (one color per roi, cycling back to the 1st if more roi
%              entries are requested than colors given -- SAME cycling convention as STYLE below, and
%              at the SAME period when left at its own default, so roi 1 and roi 5 always look
%              identical). FIT-type entries (peak contours) are NOT part of this cycle -- see STYLE
%              below for their own fixed gold/red convention. A single RGB triplet broadcasts to every
%              roi outline (today's original behavior -- style alone would then distinguish them).
%              Default a 4-color palette -- red/cyan/magenta/green -- chosen to avoid yellow (stays
%              visible across the WHOLE grayscale range a vessel patch spans, not just against its
%              darker regions -- unlike yellow, which can wash out against bright regions) and to stay
%              mutually distinguishable from one another too.
%       .exclCnsr  logical, default true -- see CENSORING below.
%       .statOverlay  struct controlling STAT-TYPE overlay rendering (see STAT-TYPE OVERLAY above) --
%              .pThresh  significance threshold on .P (default 0.05, uncorrected -- see STAT-TYPE above).
%              .colormap colormap name for .F magnitude (default 'parula').
%              .loPct/.hiPct  percentile display-range bounds on .F (default 1/99, see
%                        drawStatMapOverlay.m's own DISPLAY RANGE note).
%
% DEFAULT OVERLAY (fld char/empty form only -- see PATH SPEC above): the FIRST (up to) 3 roi labels
% already present on vessel.(fld).rois, in their own stored order -- resolved independently PER VESSEL
% when vessel is scalar, but resolved ONCE from vessel(1) and shared by every vessel in the GRID LAYOUT
% case (more than one vessel) -- see GRID LAYOUT above. A requested label not yet present is derived ON
% THE FLY via modifyRoi2.m's standard peakVox/dilate*/tissue family (idempotent -- cheap no-op if
% already present, SAME mechanism fitVessel.m/fitPatchVessels.m use), on a LOCAL COPY only -- showVessel
% never mutates or returns a vessel, so a derived label can never leak back out. This SAME on-the-fly
% derivation also applies to an EXPLICIT PATH SPEC entry naming an as-yet-undrawn roi label. A label
% that's still missing after that (not derivable) is skipped with a warning, not an error -- this is a
% display convenience, not a calculation.
%
% CENSORING (opts.exclCnsr, default true) -- frames flagged in EITHER vessel.tsCnsr.vec{r} OR
% vessel.tsCnsr_mainClust.vec{r} (per-run, true=exclude -- attachVesselMeta.m/getCnsr.m's own
% convention) are excluded from the time-/run-average this panel displays, same
% fail-open-if-absent-or-length-mismatched policy as fitIRF.m's own runCnsr (a vessel with neither
% field, e.g. a mouse vessel, is simply never censored). Set opts.exclCnsr=false to see the RAW,
% uncensored average instead (e.g. to visually confirm WHY a frame got censored).
%
% STYLE:
%   ROI-type entries (drawMaskOutline.m, cycling back to the 1st style if more than 4 are requested,
%   counted among ROI-type entries ONLY -- a FIT-type entry interleaved in the same list does not
%   consume a style/color slot):
%     1st roi : thick solid  (LineWidth 2, '-')
%     2nd roi : thick dashed (LineWidth 2, '--')
%     3rd roi : thin solid   (LineWidth 1, '-')
%     4th roi : thin dashed  (LineWidth 1, '--')
%   FIT-type entries (drawGaussianPeakContour.m): a FIXED gold/red convention, NOT cycled through
%   opts.color and NOT part of the dashed-style rotation -- SAME MAIN/SECONDARY PEAK visual vocabulary
%   fitPatchVesselsDiag.m already established (gold = the MAIN/1st peak, red = any other/SECONDARY
%   peak), so a peak contour reads consistently whether it's shown via that diagnostic figure or here.
%
% LEGEND -- exactly ONE legend PER FIGURE, positioned OUTSIDE the plot area (never overlapping a tile),
% listing the GLOBAL UNION of every overlay successfully drawn on ANY panel ACROSS THE WHOLE CALL --
% not just the panels on that particular figure (2026-08-16 -- grid overflow spans several figures, each
% with its own separate legend, see below; a label present only on vessels landing in figure 1 now still
% shows up in figure 2's own legend too, as long as figure 2 has a legend of its own at all -- Seb's own
% report: a roi ('excludeManualSinus') that happened to exist only on figure 1's vessels was silently
% missing from figure 2's legend even though figure 2 legends its own local rois just fine). Style/color
% assignment is likewise a FIXED function of a roi's own position in the requested overlay list (see
% showOneVessel's own roiOrderOf), not how many EARLIER entries happened to resolve successfully on a
% given panel -- so the same roi label is guaranteed the same color/style on every panel it's drawn on,
% never drifting out of sync with whichever panel's own handle the legend swatch happens to be built
% from. A FIT-type entry's own legend text is 'main peak (10%/90% of peak)' (its resolved index == 1) or
% 'peak N (10%/90% of peak)' (any other index) -- see resolveVesselPathSpec.m's own idx OUTPUT; the
% level-fraction text itself comes straight from drawGaussianPeakContour.m's own levelFracs OUTPUT
% (contourLevelFracStr.m), so it never drifts out of sync with what's actually drawn.
% Each figure's legend is still visually anchored to (attached via its `ax` argument to) THAT figure's
% own LAST successfully-drawn panel -- purely a positioning choice, harmless even though most of its
% entries' actual graphics handles live on a DIFFERENT figure's axes entirely (legend() labels whatever
% handles it's given, regardless of which axes they belong to -- this is what makes the GLOBAL UNION
% above possible at all: a figure with no local occurrence of some OTHER figure's roi can still legend
% it, borrowing that other figure's own handle purely for its rendering properties). Grid overflow (>15
% vessels) spans several figures -- each still gets its own single legend, now sharing identical content.
% Omitted entirely for a figure where NO panel actually resolved an overlay AT ALL (e.g. every
% requested roi label missing/unresolvable) -- a page with nothing of its own to point a legend at
% stays legend-less even if some OTHER figure has content. A mask that resolves but is null/all-false
% (nothing currently excluded/included) is NOT this case (2026-08-17) -- drawMaskOutline.m always
% returns a legend-able ghost handle for it (see that file's own OUTPUT note), so the roi still gets
% its own entry, labeled but drawing nothing on the image itself -- letting a viewer tell "this roi
% exists and is empty" apart from "this roi was never requested/found here at all".
% The axes TITLE is "<sId> <label> (<fld>)" on EVERY panel regardless, so which data field is showing is
% never ambiguous even when composing into a caller-supplied opts.ax.
%
% NO-ARG CALL -- showVessel() (zero input arguments) does not plot anything: it PRINTS a reference
% listing of every field/allowed-value (see PRINTOPTSHELP), then returns opts, fully populated with
% every default value (see DEFAULTOPTS) -- same convention as fitVessel.m's own NO-ARG CALL. Calling
% showVessel(vessel) (or any other arity) with opts omitted fills in the SAME defaults internally
% (silently -- no printing).
%
% OUTPUT
%   ax : axes handle, [1 x numel(vessel)] if vessel had more than one element (a vessel skipped for
%        missing image data gets an invalid/default graphics-object placeholder at its own index).
%        (showVessel() called with NO input arguments is the one exception -- see NO-ARG CALL above: it
%        returns opts instead.)
%
%   showVessel(vessel);                                             % default fld, first 3 rois
%   showVessel(vessel, 'tsAlignImCom');                              % explicit fld, still first 3 rois
%   showVessel(vessel, {'tsAlignImCom.rois.lumenGauss', 'tsAlignImCom.rois.lumenGauss_tissue'});
%   showVessel(vessel, {'tsIm.rois(1)', 'tsIm.patchVesselsSimultaneous.fits', 'tsIm.patchVesselsSimultaneous.fits(2)'});
%   showVessel(vessel, {'tsIm.rois(1)', 'tsIm.irf.irfIm.stats.perVoxel'});           % roi + activation map
%   ax = showVessel(vessel(3), 'tsIm', struct('ax',someAxes));       % compose into someAxes
%   ax = showVessel(vessel(3), {'tsIm.irf.irfIm.stats.perVoxel'}, ...
%           struct('ax',someAxes, 'statOverlay',struct('pThresh',0.01,'colormap','hot')));

    if nargin==0
        printOptsHelp();
        ax = defaultOpts();
        return
    end
    if nargin<2; fld = ''; end
    if nargin<3; opts = struct(); end
    opts = fillOptsDefaults(opts);

    % A dotted char (e.g. 'tsIm.rois.lumenManual') can never be a genuine vessel.(fld) field name --
    % struct field names never contain '.' -- so it can only be a PATH SPEC typed without its wrapping
    % {} -- auto-promote rather than silently resolving to a nonexistent field (see file header INPUT
    % fld note).
    if ischar(fld) && contains(fld, '.'); fld = {fld}; end

    wasCell = iscell(vessel);
    if wasCell; cellSz = size(vessel); vessel = reshape([vessel{:}], cellSz); end
    nV = numel(vessel);
    assert(isempty(opts.ax) || nV==1, 'showVessel:axNeedsScalarVessel', ...
        'opts.ax may only be given for a SCALAR vessel (got %d) -- one axes cannot host more than one vessel.', nV);

    % PATH SPEC form (fld is a cellstr) vs. today's ORIGINAL form (fld empty/char) -- see file header
    % PATH SPEC. PATH SPEC form: background field AND overlay list are BOTH one shared value for the
    % whole call (every entry's own root must match the first, asserted) -- ORIGINAL form: fld keeps
    % resolving INDEPENDENTLY per vessel (unchanged), only DEFAULT OVERLAY is forced shared in GRID
    % LAYOUT (unchanged).
    usePathSpec = iscell(fld);
    overlaySpecsShared = {}; fldShared = ''; roiShared = {};
    if usePathSpec
        overlaySpecsShared = fld;
        for i = 1:numel(overlaySpecsShared)
            [thisRoot,~,~] = resolveVesselPathSpec(vessel(1), overlaySpecsShared{i});
            if i==1; fldShared = thisRoot; end
            assert(strcmp(thisRoot,fldShared), 'showVessel:mixedFld', ...
                ['every entry of a cellstr fld argument must share the same root field as the first ' ...
                 '(''%s'') -- got ''%s'' in ''%s''.'], fldShared, thisRoot, overlaySpecsShared{i});
        end
        if isempty(overlaySpecsShared)
            fldShared = resolveDefaultFld(vessel(1));   % explicit {} -- see file header PATH SPEC
        end
    elseif nV>1
        fld1 = fld; if isempty(fld1); fld1 = resolveDefaultFld(vessel(1)); end
        roiShared = defaultRoiLabels(vessel(1), fld1);
    end

    nRow = 3; nCol = 5; nPerFig = nRow*nCol;   % see file header GRID LAYOUT
    tl = []; curFigIdx = 0;
    ax = gobjects(1,nV);
    % ONE legend per figure (see file header LEGEND note), every one of them showing the SAME GLOBAL
    % UNION of overlays drawn on ANY panel across the WHOLE call -- NOT just the panels on that
    % particular figure (2026-08-16 fix, Seb's own report: grid overflow can span several figures, and a
    % label that only happens to appear on vessels in figure 1 was previously missing from figure 2's
    % own legend entirely, even though figure 2 has a legend of its own for whatever overlays ARE local
    % to it). A legend handle is safe to reuse across figures regardless of which one its own graphics
    % object actually lives on -- legend() only copies the handle's rendering properties for the icon,
    % it does not require the handle to share the target axes (see file header LEGEND's own "harmless
    % even though some of its entries' actual graphics handles may live on an earlier panel's axes"
    % note, which already relied on exactly this). globalH/globalLbl accumulate that union (first-seen
    % order); figAnchor tracks, per group key (0 for the scalar/opts.ax cases; figIdx for the grid
    % case), ONLY the positioning anchor (.ax/.tl, updated to the LATEST successfully-drawn panel in
    % THAT figure) -- a figure with no anchor at all never drew anything and stays legend-less entirely,
    % same "omit if nothing was drawn on this figure" convention as before this fix.
    globalH = gobjects(1,0); globalLbl = {};
    figAnchor = containers.Map('KeyType','double', 'ValueType','any');
    for v = 1:nV
        if usePathSpec
            fldV = fldShared;
        else
            fldV = fld; if isempty(fldV); fldV = resolveDefaultFld(vessel(v)); end
        end
        idStr = strtrim([char(string(vessel(v).sId)) ' ' char(string(vessel(v).label))]);

        imS = imField(vessel(v), fldV);
        if isempty(imS)
            fprintf('showVessel: vessel %d/%d (%s) skipped (no %s data)\n', v, nV, idStr, fldV);
            continue
        end

        if usePathSpec
            overlaySpecsV = overlaySpecsShared;
        else
            if nV>1; roiV = roiShared; else; roiV = defaultRoiLabels(vessel(v), fldV); end
            overlaySpecsV = cellfun(@(l) [fldV '.rois.' l], roiV, 'UniformOutput',false);
        end

        groupKey = 0; tlV = [];
        if ~isempty(opts.ax)
            ax(v) = opts.ax;   % nV==1, guaranteed by the assert above
        elseif nV==1
            fig = figure; ax(v) = axes('Parent',fig);
        else
            figIdx  = ceil(v/nPerFig);
            slotIdx = mod(v-1,nPerFig) + 1;
            if figIdx ~= curFigIdx   % a pure function of v -- robust to a skipped 1st vessel of ANY batch
                fig = figure;
                tl = tiledlayout(fig, nRow, nCol, 'Padding','tight', 'TileSpacing','tight');
                curFigIdx = figIdx;
            end
            ax(v) = nexttile(tl, slotIdx);
            groupKey = figIdx;
            tlV = tl;
        end
        [legHandles, legLabels] = showOneVessel(ax(v), vessel(v), fldV, overlaySpecsV, idStr, opts.color, opts.exclCnsr, opts.statOverlay);
        if ~isempty(legHandles)
            figAnchor(groupKey) = struct('ax',ax(v), 'tl',tlV);   % anchor moves to the latest drawn panel in this figure
            for k = 1:numel(legLabels)
                if ~any(strcmp(globalLbl, legLabels{k}))   % GLOBAL union, not overwrite -- see note above
                    globalH(end+1)   = legHandles(k);   %#ok<AGROW>
                    globalLbl{end+1} = legLabels{k};   %#ok<AGROW>
                end
            end
        end
    end

    grpKeys = keys(figAnchor);
    for i = 1:numel(grpKeys)
        anc = figAnchor(grpKeys{i});
        lgd = legend(anc.ax, globalH, globalLbl, 'Interpreter','none');   % SAME global content on every figure
        if ~isempty(anc.tl)
            lgd.Layout.Tile = 'east';   % GRID LAYOUT -- spans the whole tiledlayout, outside every tile
        else
            lgd.Location = 'eastoutside';   % scalar/opts.ax -- shrinks that one axes to make room
        end
    end
end

% ---------------------------------------------------------------------------
function printOptsHelp()
    % Reference listing of every field/subfield and its allowed value(s) -- printed by the NO-ARG CALL
    % (see file header) only; a real call (any nargin>0) fills in the same defaults silently.
    fprintf('showVessel opts -- allowed value(s) per field/subfield (default in parentheses):\n');
    fprintf('  fld          : char (background field, first 3 rois auto-overlaid) | cellstr of PATH SPECs (explicit roi/fit overlay, see file header)  (last field starting with ''ts'', first 3 rois)\n');
    fprintf('  opts.ax      : axes handle (ONLY valid for a scalar vessel)         ([]; new figure+axes per vessel)\n');
    fprintf('  opts.color   : RGB triplet | Nx3 | cell of RGB triplets, cycled per ROI-type overlay (FIT-type uses a fixed gold/red convention)  (red/cyan/magenta/green)\n');
    fprintf('  opts.exclCnsr: true | false -- exclude tsCnsr/tsCnsr_mainClust frames (see CENSORING)  (true)\n');
    fprintf('  opts.statOverlay.pThresh  : significance threshold on a STAT-TYPE overlay''s own .P, uncorrected  (0.05)\n');
    fprintf('  opts.statOverlay.colormap : colormap name for a STAT-TYPE overlay''s own .F magnitude  (''parula'')\n');
    fprintf('  opts.statOverlay.loPct/.hiPct : percentile display-range bounds on .F  (1/99)\n');
end

% ---------------------------------------------------------------------------
function opts = defaultOpts()
    % Every opts field, fully populated with its default value -- see file header NO-ARG CALL/INPUT.
    opts = struct();
    opts.ax       = [];
    opts.color    = {[1 0 0], [0 0.8 0.8], [1 0 1], [0.2 0.85 0.2]};   % red / cyan / magenta / green
    opts.exclCnsr = true;
    opts.statOverlay = struct('pThresh',0.05, 'colormap','parula', 'loPct',1, 'hiPct',99);
end

% ---------------------------------------------------------------------------
function opts = fillOptsDefaults(opts)
    % Fill in any opts field left unset (missing, or [] for .color -- same "use the default" convention
    % as fitVessel.m's own fillOptsDefaults). opts.ax's own default ([]) is indistinguishable from "not
    % given" either way, so it needs no special-casing here. Generic fill-loop lives in
    % fillOptsFromDefaults.m (shared across the several files using this same pattern) -- this function
    % only layers this file's own field-specific normalization/asserts.
    opts = fillOptsFromDefaults(opts, defaultOpts());
    opts.color = normalizeColors(opts.color);   % RGB triplet | Nx3 | cell -> uniform cell of RGB triplets

    % opts.statOverlay -- nested struct, same manual per-subfield fill convention as fitIRF.m's own
    % opts.stats (fillOptsFromDefaults only fills TOP-level fields; a caller-supplied opts.statOverlay
    % missing just one subfield would otherwise leave it unset instead of defaulted).
    def = defaultOpts();
    if ~isfield(opts,'statOverlay') || isempty(opts.statOverlay); opts.statOverlay = struct(); end
    assert(isstruct(opts.statOverlay), 'showVessel:statOverlay', 'opts.statOverlay must be a struct.');
    if ~isfield(opts.statOverlay,'pThresh')  || isempty(opts.statOverlay.pThresh);  opts.statOverlay.pThresh  = def.statOverlay.pThresh;  end
    if ~isfield(opts.statOverlay,'colormap') || isempty(opts.statOverlay.colormap); opts.statOverlay.colormap = def.statOverlay.colormap; end
    if ~isfield(opts.statOverlay,'loPct')    || isempty(opts.statOverlay.loPct);    opts.statOverlay.loPct    = def.statOverlay.loPct;    end
    if ~isfield(opts.statOverlay,'hiPct')    || isempty(opts.statOverlay.hiPct);    opts.statOverlay.hiPct    = def.statOverlay.hiPct;    end
end

% ---------------------------------------------------------------------------
function colors = normalizeColors(color)
    if iscell(color)
        colors = color;
    elseif isnumeric(color) && size(color,2)==3
        colors = num2cell(color,2);   % Nx1 cell, each cell a 1x3 row -- a 1x3 input becomes a 1-element cell
    else
        error('showVessel:badColor', 'opts.color must be an RGB triplet, an Nx3 matrix, or a cell array of RGB triplets.');
    end
end

% ---------------------------------------------------------------------------
function [legHandles, legLabels] = showOneVessel(axH, vesselV, fldV, overlaySpecs, idStr, colors, exclCnsr, statOverlay)
    % Returns this panel's own legend handles/labels rather than drawing a legend itself -- the
    % caller (showVessel) draws ONE, linked to the last successfully-drawn panel per figure (see its
    % own LEGEND note) -- every panel would otherwise show an identical, redundant legend.
    legHandles = gobjects(1,0); legLabels = {};
    imS = imField(vesselV, fldV);
    voxSz = imgVoxSz(imS);
    voxSz2 = [voxSz(1) voxSz(min(2,numel(voxSz)))];   % in-plane [row col] mm -- matches fitVessel.m's own convention

    im = imS.im; if ~iscell(im); im = {im}; end
    im = censorMaskedIm(vesselV, im, exclCnsr);   % see file header CENSORING
    imAvg = mean(cat(4,im{:}), 4, 'omitnan');   % time-/run-averaged patch, matches extractMeanIm's own convention

    plotVesselPatch(axH, imAvg, voxSz2);   % see plotVesselPatch.m -- shared vessel-patch-panel primitive
    title(axH, sprintf('%s (%s)', idStr, fldV), 'Interpreter','none');

    if isempty(overlaySpecs); return; end

    % On-the-fly roi derivation (peakVox/dilate*/tissue family), LOCAL COPY only -- see file header
    % DEFAULT OVERLAY note. SAME trigger condition this file has always used (today's original
    % ~all(ismember(...)) check, not a call on every '.rois' spec regardless of need -- that family
    % (esp. 'tissue', a full getRoiBckgrndMask.m fit) is NOT free, so this stays conditional): only when
    % at least one requested, already-LABELED ('<fldV>.rois.<label>' form) roi isn't yet present. A
    % '<fldV>.rois(N)' numeric-index spec never triggers this -- nothing to derive by position.
    roiLabelPrefix = [fldV '.rois.'];
    roiLabelSpecs = overlaySpecs(startsWith(overlaySpecs, roiLabelPrefix));
    if ~isempty(roiLabelSpecs)
        requestedLabels = cellfun(@(s) s(numel(roiLabelPrefix)+1:end), roiLabelSpecs, 'UniformOutput',false);
        haveLabels = {};
        if isfield(vesselV.(fldV),'rois') && ~isempty(vesselV.(fldV).rois)
            haveLabels = {vesselV.(fldV).rois.label};
        end
        if ~all(ismember(requestedLabels, haveLabels))
            vesselV = modifyRoi2(vesselV, fldV, {'peakVox','dilate1','dilate1p5','dilate2','tissue'});
        end
    end

    goldCol = [1 0.85 0.1]; redCol = [1 0 0];   % SAME convention as fitPatchVesselsDiag.m's own
                                                 % MAIN/SECONDARY PEAK colors -- see file header STYLE.
    styles = {{2,'-'}, {2,'--'}, {1,'-'}, {1,'--'}};   % ROI-type entries only -- thick solid / thick
                                                        % dashed / thin solid / thin dashed
    % roiOrderOf -- each spec's own FIXED position among ROI-type entries in overlaySpecs itself
    % (1-based; 0 for a non-roi spec), computed ONCE from the spec strings before any resolving happens
    % -- see file header STYLE's own "cycling ... counted among ROI-type entries ONLY". Deliberately NOT
    % a running counter incremented only on successful resolution (2026-08-16 fix, Seb's own report: a
    % roi missing on THIS panel but present on others used to shift every SUBSEQUENT roi's color/style
    % out of sync between panels -- e.g. panel A has rois {lumen, exclude} -> exclude gets style/color
    % #2, but panel B is missing lumen -> exclude becomes THIS panel's #1 -> style/color #1 instead. The
    % shared legend swatch (built from whichever panel's handle got captured first) then visibly
    % mismatched whichever OTHER panel's own outline actually used the other color). A roi-type spec
    % names a '.rois' segment (see resolveVesselPathSpec.m's own SYNTAX, both '.rois.<label>' and
    % '.rois(N)' forms) -- knowable from the spec STRING alone, no need to resolve it first.
    isRoiSpec = startsWith(overlaySpecs, [fldV '.rois']);
    roiOrderOf = cumsum(isRoiSpec) .* isRoiSpec;
    [ny,nx] = size(imAvg);
    XIsoFit = []; YIsoFit = [];   % built lazily, ONCE, only if a FIT-type entry is actually requested --
                                   % see drawGaussianPeakContour.m's own "caller builds the grid once" note.

    for k = 1:numel(overlaySpecs)
        spec = overlaySpecs{k};
        try
            [thisRoot, node, idx] = resolveVesselPathSpec(vesselV, spec);
        catch ME
            warning('showVessel:overlaySpecNotFound', ...
                'vessel %s: overlay spec ''%s'' could not be resolved (%s) -- skipped.', idStr, spec, ME.message);
            continue
        end
        assert(strcmp(thisRoot,fldV), 'showVessel:specRootMismatch', ...
            'overlay spec ''%s'' does not share this panel''s own background field ''%s''.', spec, fldV);

        if isstruct(node) && isfield(node,'mask')
            % ROI-type -- st/col from this spec's own FIXED roiOrderOf(k) (precomputed above from the
            % spec LIST itself), not a per-panel running counter -- see that comment for why.
            st  = styles{mod(roiOrderOf(k)-1,numel(styles))+1};
            col = colors{mod(roiOrderOf(k)-1,numel(colors))+1};
            % drawMaskOutline.m always returns a valid (possibly ghost, for a null/all-false mask)
            % handle now (2026-08-17) -- so a roi with nothing currently excluded still gets its own
            % legend entry, rather than being silently omitted as if it didn't exist at all.
            h = drawMaskOutline(axH, node.mask, voxSz2, col, st{2}, st{1});
            legHandles(end+1) = h;   %#ok<AGROW>
            legLabels{end+1}  = node.label;   %#ok<AGROW>
        elseif isstruct(node) && isfield(node,'a') && isfield(node,'x0') && (isfield(node,'sx') || isfield(node,'radius'))
            % FIT-type (one fitPatchVessels.m peak, gaussianModel.m's own param convention) -- see file
            % header PATH SPEC/STYLE. Checks BOTH shape-field conventions (2026-08-18) -- this project's
            % own fitVessel.m/fitPatchVessels.m ALWAYS use radius/aspectRatio, never .sx/.sy (see
            % drawGaussianPeakContour.m's own header, which already handles that conversion correctly);
            % this detection check had simply never been updated to match, so a real fit-peak overlay
            % silently fell through to the "not drawable" warning below instead of ever reaching that
            % already-correct drawing code.
            if isempty(XIsoFit)
                [XIsoFit,YIsoFit] = buildIsochromatGrid(7, nx, ny, (nx+1)/2, (ny+1)/2, voxSz2);
            end
            col = redCol; if idx==1; col = goldCol; end
            [h,~,levelFracs] = drawGaussianPeakContour(axH, node, XIsoFit, YIsoFit, col, 0.5);
            lbl = 'main peak'; if idx~=1; lbl = sprintf('peak %d', idx); end
            lbl = sprintf('%s (%s of peak)', lbl, contourLevelFracStr(levelFracs));
            legHandles(end+1) = h;   %#ok<AGROW>
            legLabels{end+1}  = lbl;   %#ok<AGROW>
        elseif isstruct(node) && isfield(node,'F') && isfield(node,'P')
            % STAT-TYPE (fitIRF.m's own per-voxel test map, e.g. '<fld>.irf.irfIm.stats.perVoxel') -- see
            % file header STAT-TYPE OVERLAY/drawStatMapOverlay.m for the drawing recipe.
            h = drawStatMapOverlay(axH, node, voxSz2, statOverlay.pThresh, statOverlay.colormap, ...
                statOverlay.loPct, statOverlay.hiPct);
            legHandles(end+1) = h;   %#ok<AGROW>
            legLabels{end+1}  = sprintf('activation map (F, P<%.3g)', statOverlay.pThresh);   %#ok<AGROW>
        else
            warning('showVessel:overlayNotDrawable', ...
                'vessel %s: overlay spec ''%s'' resolved to something that isn''t a drawable roi or fit -- skipped.', idStr, spec);
        end
    end
end

% ---------------------------------------------------------------------------
function labels = defaultRoiLabels(vesselV, fld)
    % Default overlay (file header DEFAULT OVERLAY): the FIRST (up to) 3 roi labels already present, in
    % stored order.
    labels = {};
    if isfield(vesselV,fld) && isstruct(vesselV.(fld)) && isfield(vesselV.(fld),'rois') && ~isempty(vesselV.(fld).rois)
        allLabels = {vesselV.(fld).rois.label};
        labels = allLabels(1:min(3,numel(allLabels)));
    end
end

% ---------------------------------------------------------------------------
function fld = resolveDefaultFld(vesselV)
    % Same default as fitVessel.m's own INPUT fld (duplicated locally -- see this codebase's
    % established small-helper-duplication convention, e.g. fitVessel.m's own gaussianDefaults note):
    % the LAST field on this vessel whose name starts with 'ts' AND actually carries image data (via
    % imField -- excludes e.g. tsCnsr_mainClust, a censor mask field with no .im at all), in
    % fieldnames(vesselV) order.
    names = fieldnames(vesselV);
    tsNames = names(startsWith(names,'ts'));
    tsNames = tsNames(~cellfun(@(n) isempty(imField(vesselV,n)), tsNames));
    assert(~isempty(tsNames), 'showVessel:noTsFld', ...
        'vessel %s: no field starting with ''ts'' with image data found -- specify fld explicitly.', char(string(vesselV.label)));
    fld = tsNames{end};
end

% ---------------------------------------------------------------------------
function S = imField(vessel, fld)
    % Duplicated from fitVessel.m -- mirrors getCom.m's own imField exactly (top-level vessel.(fld)
    % preferred, legacy nested vessel.im.(fld) fallback).
    S = [];
    if isfield(vessel,fld) && isstruct(vessel.(fld)) && isfield(vessel.(fld),'im') && ~isempty(vessel.(fld).im)
        S = vessel.(fld);
    elseif isfield(vessel,'im') && isstruct(vessel.im) && isfield(vessel.im,fld) && ...
           isstruct(vessel.im.(fld)) && isfield(vessel.im.(fld),'im') && ~isempty(vessel.im.(fld).im)
        S = vessel.im.(fld);
    end
end

% ---------------------------------------------------------------------------
function voxSz = imgVoxSz(imStruct)
    % Duplicated from fitVessel.m.
    assert(isfield(imStruct,'vSize') && ~isempty(imStruct.vSize), 'showVessel:noVoxSize', ...
        'image struct is missing .vSize (mm/voxel); attach it before plotting.');
    voxSz = imStruct.vSize(:);
end

% ---------------------------------------------------------------------------
function im = censorMaskedIm(vesselV, im, exclCnsr)
    % NaN-out censored frames (dim 4 of each per-run array, BEFORE cat(4,...)+mean(...,'omitnan')) --
    % see file header CENSORING. Excludes a frame if EITHER vessel.tsCnsr.vec{r} OR
    % vessel.tsCnsr_mainClust.vec{r} flags it (true=exclude -- attachVesselMeta.m/getCnsr.m's own
    % convention). Fail-open (run left untouched) if a censor field is absent or its length doesn't
    % match that run's own frame count -- same convention as fitIRF.m's own runCnsr.
    if ~exclCnsr; return; end
    cnsrA = cnsrVecLocal(vesselV, 'tsCnsr');
    cnsrB = cnsrVecLocal(vesselV, 'tsCnsr_mainClust');
    if isempty(cnsrA) && isempty(cnsrB); return; end
    for r = 1:numel(im)
        T = size(im{r}, 4);
        keep = runKeepLocal(cnsrA, r, T) & runKeepLocal(cnsrB, r, T);
        if ~all(keep)
            im{r}(:,:,:,~keep) = NaN;
        end
    end
end

% ---------------------------------------------------------------------------
function cnsr = cnsrVecLocal(vesselV, fldName)
    % Duplicated from fitIRF.m's own vesselCnsr -- {} if the field is absent (e.g. a mouse vessel).
    cnsr = {};
    if isfield(vesselV,fldName) && isstruct(vesselV.(fldName)) && isfield(vesselV.(fldName),'vec')
        cnsr = vesselV.(fldName).vec;
        if ~iscell(cnsr); cnsr = {cnsr}; end
    end
end

% ---------------------------------------------------------------------------
function keep = runKeepLocal(cnsr, r, T)
    % Duplicated from fitIRF.m's own runCnsr -- [1 x T] KEEP vector (inverted from cnsr's own
    % true=exclude convention); all-true (nothing excluded) if this run's own censor vector is
    % missing or its length doesn't match T -- fail open, not an error.
    keep = true(1,T);
    if isempty(cnsr) || numel(cnsr) < r || isempty(cnsr{r}) || numel(cnsr{r}) ~= T; return; end
    keep = ~logical(cnsr{r}(:))';
end
