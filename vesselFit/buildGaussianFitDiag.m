function [diagS, titleStr] = buildGaussianFitDiag(vessel, fld, outFld, run, frame)
% BUILDGAUSSIANFITDIAG  Reconstruct the DIAG STRUCT fitVessel.m's own opts.verbose would have built
% internally (see fitVesselDiag.m's own file header DIAG STRUCT) from a stored vessel.(fld).(outFld)
% result -- recomputing validMask/predicted IDENTICALLY from the stored fit parameters + mask
% provenance, never from a re-run of the optimizer.
%
% Promoted out of fitVesselDiag.m's own local buildDiagFromVessel (2026-08-18), unchanged in substance,
% so a SECOND renderer -- plotGaussianFitPanels.m, a curated manuscript-example figure, not a QC tool --
% can reuse the exact same reconstruction rather than re-deriving/duplicating it (this file family has
% already been bitten, repeatedly, by hand-reconstructed copies of fit-provenance math drifting stale --
% see [[tool_vfmritools_vessel_figure_terminology]]'s own "KEEP IN SYNC" history). fitVesselDiag.m now
% calls this function instead of carrying its own local copy.
%
% See buildGaussianFitDiagFromPeak.m for the SIBLING constructor -- same DIAG STRUCT shape, built from
% ONE peak of a fitPatchVessels.m result (vessel.(fld).(outFld).fits(k)) instead of a flat
% fitVessel.m-shaped result.
%
% MULTI-COMPONENT RESULTS (2026-09-04) -- also accepts a fitVesselPatchTimeSeries.m result, whose fit
% is a SUM of components (K gaussian peaks plus one background). Such a result carries a flat
% main-peak alias, so every field below that describes "the fit" keeps its old single-peak meaning and
% every existing consumer is unaffected. THREE extra fields carry what the flat view structurally
% cannot:
%   .nComponent   how many components the fit actually had (1 for a fitVessel.m-shaped result).
%   .modelFull    [ny x nx] -- EVERY component summed. This is the honest predicted image. Rebuilding
%                 it from .fitVals alone would give a MAIN-PEAK-ONLY model and therefore render every
%                 secondary peak as residual, which is the specific defect these fields exist to fix.
%   .otherFull    [ny x nx] -- the SECONDARY peaks alone, background EXCLUDED. Subtract it from the
%                 measured patch to get the main-peak-only signal a 1-D radial profile is about (the
%                 same "secondary subtracted" treatment fitPatchVesselsDiag.m already applies to its
%                 own profile/scatter panels).
%   .otherFitVals 1xM struct array, one per secondary peak, for drawing their contours. NOTE these
%                 carry the REAL background in .b (a contour level is b + a*frac) even though the
%                 .otherFull SUM excludes it -- two different uses needing two different b.
% .predicted (masked) is likewise the FULL model, so it means the same thing for both kinds of result.
%
%   vessel : ONE scalar vessel struct -- see fitVesselDiag.m's own vessel-mode INPUT note.
%   fld    : the field fitVessel.m was called against (e.g. 'tsAlignImAfni').
%   outFld : fitVessel.m's own opts.outFld (e.g. 'gaussFitXY') -- vessel.(fld).(outFld) must exist.
%   run    : which run to diagnose -- see fitVesselDiag.m's own opts.run note.
%   frame  : which frame within that run -- see fitVesselDiag.m's own opts.frame note.
%
%   [diagS, titleStr] = buildGaussianFitDiag(vesselHuman{3}, 'tsAlignImAfni', 'gaussFitXY', 1, 5);

    % Checked FIRST, separately from the isfield(vessel,fld) check below: isfield() on a non-struct
    % (e.g. the whole vesselHuman CELL ARRAY, passed instead of one indexed vesselHuman{k}) silently
    % returns false rather than erroring -- without this check, that mistake surfaces as the field-check
    % assert below instead, misleadingly blaming "vessel.<fld> is missing" when the real problem is
    % vessel itself being the wrong shape (see [[feedback_no_fallbacks]] -- a confusing error is exactly
    % what that pattern exists to prevent).
    assert(isstruct(vessel) && isscalar(vessel), 'buildGaussianFitDiag:notScalarVessel', ...
        ['vessel must be ONE scalar vessel struct, not an array or cell -- index into it yourself ' ...
         '(e.g. vesselHuman{1} or vesselHuman(1)).']);
    assert(isfield(vessel,fld) && isstruct(vessel.(fld)) && isfield(vessel.(fld),'im') && ~isempty(vessel.(fld).im), ...
        'buildGaussianFitDiag:noFld', 'vessel.%s is missing or has no .im.', fld);
    imS = vessel.(fld);
    assert(isfield(imS,outFld) && isstruct(imS.(outFld)), 'buildGaussianFitDiag:noOutFld', ...
        'vessel.%s.%s is missing -- run fitVessel(vessel, ''%s'', ''gaussian''|''parabolic'', struct(''outFld'',''%s'',...)) first.', ...
        fld, outFld, fld, outFld);
    res = imS.(outFld);
    assert(isfield(res,'method') && any(strcmp(res.method,{'gaussian','parabolic'})), 'buildGaussianFitDiag:notModelFit', ...
        'vessel.%s.%s.method is ''%s'' -- must be ''gaussian'' or ''parabolic''.', fld, outFld, char(string(res.method)));

    gspec = methodSpecMinimal(res.method);

    im = imS.im; wasImCell = iscell(im); if ~wasImCell; im = {im}; end
    nRunRaw = numel(im);

    % Does "run" still mean the same thing on the raw image as it did at fit time? packFitOut
    % cell-wraps every per-parameter field (one cell per pseudo-run) only when the fit itself saw >1
    % pseudo-run -- see fitVessel.m's own packFitOut. If that pseudo-run count doesn't match
    % vessel.(fld).im's own raw run count, opts.runAvg/opts.runCat must have combined runs into
    % something other than "one pseudo-run per original run", and there is no stored provenance to
    % reconstruct exactly what that combination was -- error loudly rather than silently draw the
    % wrong picture (see [[feedback_no_fallbacks]]).
    % Any parameter that varies between runs is a {1 x nRun} cell; one identical across runs is
    % stored bare (fitVesselPatchTimeSeries.m's collapseIdentical, 2026-09-10), so the fit's run count
    % is the widest cell over ALL parameters, and a result with no cell at all (every parameter
    % per-vessel) carries no run count of its own to disagree with the raw image.
    nRunAtFit = 1;
    for p = 1:numel(gspec.paramNames)
        v = res.(gspec.paramNames{p});
        if iscell(v); nRunAtFit = max(nRunAtFit, numel(v)); end
    end
    assert(nRunAtFit == 1 || nRunAtFit == nRunRaw, 'buildGaussianFitDiag:runMismatch', ...
        ['vessel.%s.im has %d run(s) but vessel.%s.%s was fit against %d pseudo-run(s) -- this ' ...
         'reconstruction assumes opts.runAvg=false/opts.runCat=false at fit time (pseudo-run == ' ...
         'original run, the only case a raw run index still means anything); a different run ' ...
         'combination was used here and can''t be reconstructed without re-running the fit.'], ...
        fld, nRunRaw, fld, outFld, nRunAtFit);
    assert(run>=1 && run<=nRunRaw, 'buildGaussianFitDiag:badRun', ...
        'run %d out of range (vessel.%s.im has %d run(s)).', run, fld, nRunRaw);

    im4 = double(im{run}(:,:,1,:));
    [ny,nx,~,Traw] = size(im4);

    if res.timeAvg
        assert(frame==1, 'buildGaussianFitDiag:timeAvgFrame', ...
            ['vessel.%s.%s was fit with opts.timeAvg=true (one time-averaged frame per run) -- frame ' ...
             'must be 1 (got %d).'], fld, outFld, frame);
        frameIm = mean(im4,4,'omitnan');
    else
        assert(frame>=1 && frame<=Traw, 'buildGaussianFitDiag:badFrame', ...
            'frame %d out of range (run %d has %d frame(s)).', frame, run, Traw);
        frameIm = im4(:,:,1,frame);
    end
    nanM = isnan(frameIm); if any(nanM(:)); frameIm(nanM) = mean(frameIm(~nanM),'omitnan'); end

    voxSz = imS.vSize(:);
    voxSz2 = [voxSz(1) voxSz(min(2,numel(voxSz)))];   % in-plane [row col] mm -- mirrors fitVessel.m's own fitOneRun
    ctrCol = (nx+1)/2; ctrRow = (ny+1)/2;
    largeSz = [ny nx];

    % Every mask label (outFld's own .includeMask/.excludeMask provenance) must already be present on
    % vessel.(fld).rois -- fitVessel.m appends any it derived on the fly and returns the vessel WITH
    % them, so a missing label here means a stale/earlier vessel snapshot was passed instead, not a
    % genuinely absent mask (unionMasksLocal/resolveIncludeMaskLocal would otherwise silently treat a
    % missing label as contributing nothing, drawing a wrong validMask with no warning at all).
    % Guarded on the label list being NON-EMPTY, and the label list read defensively (2026-09-02):
    % a fit that used NO mask at all has nothing to check, and requiring .rois anyway made an
    % unmasked fit undiagnosable. Building ismember's own second argument dereferenced imS.rois
    % EAGERLY, so a fld with no .rois field threw "Unrecognized field name rois" right here even
    % though both label lists were {} -- hit for real on vessel.tsRawIm (nothing is ever drawn on the
    % raw patches). Same fix, same reason, as fitVessel.m's own wantLabels guard; and the mask
    % builders just below (resolveIncludeMaskLocal/unionMasksLocal) already short-circuit on an empty
    % list before touching .rois, so this assert was the only thing standing in the way.
    wantMasks = [res.includeMask, res.excludeMask];
    if ~isempty(wantMasks)
        haveLabels = {};
        if isfield(imS,'rois') && ~isempty(imS.rois); haveLabels = {imS.rois.label}; end
        assert(all(ismember(wantMasks, haveLabels)), 'buildGaussianFitDiag:missingRois', ...
            ['vessel.%s.rois is missing one or more of the mask label(s) vessel.%s.%s was actually fit ' ...
             'with (.includeMask/.excludeMask provenance) -- pass the SAME vessel fitVessel.m returned ' ...
             '(it appends any on-the-fly-derived labels), not an earlier snapshot.'], fld, fld, outFld);
    end
    includeMaskLarge = resolveIncludeMaskLocal(vessel, fld, res.includeMask, largeSz);
    excludeMaskLarge = unionMasksLocal(vessel, fld, res.excludeMask, largeSz);
    validMask = includeMaskLarge(:,:,1) & ~excludeMaskLarge(:,:,1);

    nIso = res.nIsochromatPerVoxDim;
    [XIso, YIso] = buildIsochromatGrid(nIso, nx, ny, ctrCol, ctrRow, voxSz2);

    nP = numel(gspec.paramNames);
    fitVals = struct();
    for p = 1:nP
        name = gspec.paramNames{p};
        v = perRunVal(res.(name), run);
        if isscalar(v); fitVals.(name) = v; else; fitVals.(name) = v(frame); end
    end

    % See toModelArgs for the radius/aspectRatio -> literal sx,sy conversion and why diagS.fitVals is
    % deliberately left UNCONVERTED.
    valsCellForModel = toModelArgs(fitVals, gspec.paramNames, res.method);

    % MAIN component, at FULL [ny nx] resolution (not masked -- a renderer needs the whole patch).
    mainFull = reduceIsochromatFull(gspec.modelFun(XIso, YIso, valsCellForModel{:}), nIso, ny, nx);

    % MULTI-COMPONENT RECONSTRUCTION (2026-09-04). Everything above reads the flat, main-peak view --
    % which a fitVesselPatchTimeSeries.m result also carries (its own flat alias), so it needs no
    % special casing at all. What that view CANNOT give is the full prediction when the fit had more
    % than one component: rebuilding from it alone yields a MAIN-PEAK-ONLY model, which then renders
    % every secondary peak as residual -- the exact gap these fields close.
    %
    % A single-component (fitVessel.m-shaped) result gets the same fields with nComponent=1, otherFull
    % all zeros and otherFitVals empty -- so a consumer never has to ask which kind of result it is
    % holding, and its own numbers are unchanged.
    otherFull = zeros(ny, nx);
    otherFitVals = struct([]);
    nComponent = 1;
    if isfield(res,'model') && isstruct(res.model)
        [otherFull, otherFitVals, nComponent] = secondaryComponents(res, run, frame, fitVals.b, ...
                                                                    XIso, YIso, nIso, ny, nx);
    end
    modelFull = mainFull + otherFull;

    diagS = struct();
    diagS.im              = frameIm;
    diagS.voxSz2          = voxSz2;
    diagS.modelFun        = gspec.modelFun;
    diagS.paramNames      = gspec.paramNames;
    diagS.shapeParamNames = gspec.shapeParamNames;
    diagS.fitVals         = fitVals;
    diagS.validMask        = validMask;
    diagS.predicted        = modelFull(validMask);   % the FULL model, masked -- see above
    diagS.nIso              = res.nIsochromatPerVoxDim;
    diagS.method            = res.method;
    % --- multi-component additions ---
    diagS.nComponent        = nComponent;    % 1 for a single-model fit
    diagS.modelFull         = modelFull;     % [ny nx], EVERY component summed (+ background)
    diagS.otherFull         = otherFull;     % [ny nx], SECONDARY peaks only, background EXCLUDED --
                                              % subtract from the measured patch to get the
                                              % main-peak-only signal the 1-D profile is about
    diagS.otherFitVals      = otherFitVals;  % 1xM struct array, one per secondary (for red contours)

    idStr = strtrim([char(string(vessel.sId)) ' ' char(string(vessel.label))]);
    titleStr = sprintf('%s (%s.%s, %s, run %d/%d, frame %d/%d)', idStr, fld, outFld, res.method, run, nRunRaw, frame, Traw);
end

% ---------------------------------------------------------------------------
function gspec = methodSpecMinimal(method)
    % paramNames/modelFun/shapeParamNames only -- mirrors fitVessel.m's own gaussianMethodSpec/
    % parabolicMethodSpec (duplicated locally rather than cross-file-called -- this pipeline's
    % established small-helper-duplication convention, see fitIRF.m's gaussianDefaults note). No
    % defaultsFun here: a post-hoc diagnostic never needs a StartPoint heuristic.
    %
    % 'gaussian's own paramNames/shapeParamNames are radius/aspectRatio (this project's fitVessel.m
    % always optimizes that pair, see its own PARAMETERIZATION note -- there is no .sx/.sy field
    % anywhere in this project) -- used to read res.(name) correctly, and to know which two fields are
    % the shape pair when a renderer needs to convert to literal sx,sy itself (contour/radial-profile
    % math) -- see that renderer's own note.
    switch method
        case 'gaussian'
            gspec.paramNames      = {'a','x0','y0','radius','aspectRatio','theta','b'};
            gspec.modelFun        = @gaussianModel;
            gspec.shapeParamNames = {'radius','aspectRatio'};
        case 'parabolic'
            gspec.paramNames      = {'a','x0','y0','Rx','Ry','theta','b','bTissue'};
            gspec.modelFun        = @parabolicModel;
            gspec.shapeParamNames = {'Rx','Ry'};
        otherwise
            error('buildGaussianFitDiag:unknownMethod', 'unknown method ''%s''.', method);
    end
end

% ---------------------------------------------------------------------------
function mUnion = unionMasksLocal(vessel, fld, labelList, largeSz)
    % Verbatim copy of fitVessel.m's own unionMasks (duplicated locally, same convention as
    % methodSpecMinimal above) -- OR-union of every ROI in labelList that actually exists on
    % vessel.(fld).rois.
    mUnion = false(largeSz);
    if isempty(labelList); return; end
    labels = {vessel.(fld).rois.label};
    for i = 1:numel(labelList)
        idx = find(strcmp(labels, labelList{i}), 1);
        if isempty(idx); continue; end
        mUnion = mUnion | embedIfSmallLocal(vessel.(fld).rois(idx).mask, vessel, largeSz);
    end
end

% ---------------------------------------------------------------------------
function m = resolveIncludeMaskLocal(vessel, fld, includeMaskArg, largeSz)
    % Verbatim copy of fitVessel.m's own resolveIncludeMask.
    if isempty(includeMaskArg); m = true(largeSz); return; end
    m = unionMasksLocal(vessel, fld, includeMaskArg, largeSz);
end

% ---------------------------------------------------------------------------
function mask = embedIfSmallLocal(mask, vessel, largeSz)
    % Verbatim copy of fitVessel.m's own embedIfSmall.
    if size(mask,1)==largeSz(1) && size(mask,2)==largeSz(2); return; end
    assert(size(mask,1)==size(mask,2), 'buildGaussianFitDiag:badMaskShape', ...
           'mask must be square (cropSz+1) to infer its final-crop-window placement.');
    cropSz = size(mask,1) - 1;
    [offRow, offCol, nR] = finalCropWindow(vessel.com, vessel.cropXlim, vessel.cropYlim, cropSz);
    assert(offRow>=0 && offCol>=0 && offRow+nR<=largeSz(1) && offCol+nR<=largeSz(2), ...
           'buildGaussianFitDiag:badCropWindow', 'final-crop window falls outside the current large image frame.');
    big = false(largeSz);
    big(offRow+1:offRow+nR, offCol+1:offCol+nR, :) = mask;
    mask = big;
end

% ---------------------------------------------------------------------------
function mv = toModelArgs(fitVals, paramNames, method)
    % A parameter struct -> the model formula's own POSITIONAL argument list. For 'gaussian', vals{4:5}
    % are radius/aspectRatio (this project's sole parameterization -- see fitVessel.m's own
    % PARAMETERIZATION note; there is no .sx/.sy field anywhere in this project), which gaussianModel.m
    % has no concept of, so they are converted to literal sx,sy here, as plain UNNAMED locals, purely to
    % evaluate the model. diagS.fitVals stays in radius/aspectRatio, UNCONVERTED -- every renderer does
    % this same conversion for its own contour/radial-profile math, so this struct is never where sx/sy
    % gets stored.
    vals = cellfun(@(n) fitVals.(n), paramNames, 'UniformOutput',false);
    if strcmp(method,'gaussian')
        sqrtAR = sqrt(vals{5});
        mv = [vals(1:3), {vals{4}*sqrtAR, vals{4}/sqrtAR}, vals(6:7)];
    else
        mv = vals;
    end
end

% ---------------------------------------------------------------------------
function [otherFull, otherFitVals, nComponent] = secondaryComponents(res, run, frame, bMain, ...
                                                                     XIso, YIso, nIso, ny, nx)
    % Every gaussian component EXCEPT the main one (array position 1), summed at full [ny nx]
    % resolution.
    %
    % BACKGROUND IS EXCLUDED FROM THE SUM, deliberately: the main-peak evaluation already carries b
    % (it arrives through the flat alias), so adding it again per secondary would multiply-count it.
    % A consumer subtracting otherFull from the measured patch also wants the secondaries gone WITHOUT
    % losing the baseline the main peak sits on. The returned otherFitVals nevertheless carry the REAL
    % background in .b, because that is a CONTOUR-drawing input (drawGaussianPeakContour.m places its
    % levels at b + a*frac) -- so the two uses need different b, and conflating them would either
    % double-count the baseline in the image or draw every secondary contour at the wrong level.
    otherFull = zeros(ny, nx);
    otherFitVals = struct([]);
    nComponent = 0;
    if ~isfield(res.model,'gaussian'); return; end
    g = res.model.gaussian;
    nComponent = numel(g);
    if isfield(res.model,'background'); nComponent = nComponent + numel(res.model.background); end
    if numel(g) < 2; return; end

    pnames = {'a','x0','y0','radius','aspectRatio','theta'};
    for k = 2:numel(g)
        fv = struct();
        for p = 1:numel(pnames)
            % {1 x nRun} cell when the parameter varies between runs, bare otherwise (perRunVal.m
            % covers both); the run-count agreement with the raw image was already asserted above.
            % Within a run a parameter is either one value (any coarser granularity) or one per frame.
            v = perRunVal(g(k).params.(pnames{p}), run);
            if isscalar(v); fv.(pnames{p}) = v; else; fv.(pnames{p}) = v(frame); end
        end
        fv.b = bMain;
        mvArgs = toModelArgs(fv, [pnames {'b'}], 'gaussian');
        mvArgs{7} = 0;   % peak only -- see BACKGROUND IS EXCLUDED above
        otherFull = otherFull + reduceIsochromatFull(gaussianModel(XIso, YIso, mvArgs{:}), nIso, ny, nx);
        if isempty(otherFitVals); otherFitVals = fv; else; otherFitVals(end+1) = fv; end %#ok<AGROW>
    end
end
