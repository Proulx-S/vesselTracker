function out = fitComponentFit(im4, validMask, voxSz2, components, opts)
% FITCOMPONENTFIT  Joint fit of a SUM OF MODEL COMPONENTS to one image patch's own frames, with
% per-parameter granularity (perVessel/perFrame/perRun/perKnot) on every parameter of every component.
% One lsqnonlin call over the whole thing.
%
%     pred(t) = sum_c  modelFun_c(X, Y, toModelParams_c(vals_c(t)){:})
%
% The shared solver behind fitVesselPatchTimeSeries.m. Generalizes fitVessel.m's own fitJointPath from
% ONE model's parameters to a LIST of components' parameters -- the free vector goes from being packed
% by (parameter, granularity-group) to (COMPONENT, parameter, granularity-group). That extra index is
% the entire difference; the granularity machinery, the JacobPattern reasoning, and the bound/start
% handling are all fitVessel.m's, unchanged in substance.
%
% A PURE NUMERICAL ENGINE -- deliberately knows nothing about vessel structs, ROI label masks, dot-path
% sources, seeding heuristics or default bounds. Every start value and bound arrives as a LITERAL
% NUMBER that the caller has already resolved. That is what makes it testable in isolation and reusable
% by callers with completely different seeding conventions; it is also why there is no .source field
% anywhere here (see fitVesselPatchTimeSeries.m, which owns seeding).
%
% EVERY FREE PARAMETER IS "nGroups COLUMNS, COLUMN g TOUCHING THE FRAMES WHERE groupIdx==g".
% This is the one idea the whole file rests on. A parameter's GRANULARITY is fully described by a
% [1 x T] frame->group vector plus that vector's own group count:
%     'perVessel' -> groupIdx = ones(1,T),  nGroups = 1        (one shared value; a dense column)
%     'perFrame'  -> groupIdx = 1:T,        nGroups = T        (independent per frame)
%     'perRun'    -> groupIdx = opts.runIdx,  nGroups = #runs
%     'perKnot'   -> groupIdx = opts.knotIdx, nGroups = #knots
% Once granularity is represented that way, packing, unpacking, bound filling and the Jacobian sparsity
% pattern are all ONE loop over slots -- fitVessel.m needs four near-identical arms for these (one per
% mode) precisely because it hardcodes the four cases instead of the vector behind them. Adding a fifth
% granularity here would mean supplying another frame->group vector and nothing else -- true for a fifth
% GROUPING granularity; mode='perRunPoly' below is different in kind (a design MATRIX, not a group
% selector) and gets its own value-computation path for exactly that reason -- see that section.
%
% mode='perRunPoly([<degrees>])' (2026-09-07) -- an INDEPENDENT low-order Legendre polynomial in
% within-run frame time, with its own coefficients PER ORIGINAL RUN, filling the gap between 'perRun'
% (one constant per run) and 'perFrame' (fully free per frame): e.g. a vessel position that drifts
% smoothly within a run, with its own drift trajectory per run. The degree list is EMBEDDED in the mode
% string itself (there is no separate opts field for it), accepted in any of these spellings --
%     'perRunPoly([1 2])'   'perRunPoly(1,3)'   'perRunPoly[1 2]'   'perRunPoly(0:4)'
% -- a comma/space-separated list of non-negative integers, each optionally a colon RANGE ('a:b' or
% 'a:step:b'), itself optionally wrapped in one more layer of ()/[] (so 'perRunPoly([1 2])' and
% 'perRunPoly(1,2)' are the same request). Degree 0 = constant/intercept, degree 1 = linear drift, etc.
% -- omitting 0 (e.g. 'perRunPoly([1 2])') is allowed and means "no intercept term at all", not "zero
% intercept" (see SEEDING below). Parsed and validated once per (component,parameter), then normalized
% to a CANONICAL spelling -- 'perRunPoly([<sorted unique degrees, space-separated>])' -- which is what
% .mode echoes back and what a caller sees on output; two textually different requests for the same
% degree set therefore compare equal once resolved (though see fitVesselPatchTimeSeries.m's own replicate
% cache key, which is generous and keys off the CALLER's original spelling instead).
%
% ONLY ESTIMATED (fixed=false) is supported -- fixed=true errors (fitComponentFit:
% perRunPolyFixedUnsupported); pin a run's own trend by fitting/precomputing it separately instead
% (mode='perRun', fixed=true, literal per-run .start).
%
% SEEDING (2026-09-07, Seb's own call): the caller still supplies exactly ONE scalar .start and ONE
% scalar .lower/.upper for the parameter, precisely as for every other mode (see PURE NUMERICAL ENGINE --
% this engine still takes no default-heuristic values of its own). Internally that scalar seeds the
% DEGREE-0 coefficient of EVERY run identically; every OTHER coefficient starts at 0 (no assumed trend).
% Bounds mirror this: the degree-0 coefficient gets the caller's own [.lower .upper] on every run; every
% other coefficient gets a DEFAULT symmetric bound of +-(.upper-.lower)/2 (letting one polynomial term
% swing across the parameter's own full caller-given dynamic range, centered on 0) -- there is no separate
% per-degree bound override channel. A degree list that omits 0 entirely gets no degree-0 seeding at all
% (every present coefficient starts at 0, symmetric bound); the caller's .start value is then unused by
% this parameter (still required, for shape/API uniformity with every other mode).
%
% OUTPUT, DELIBERATELY DIFFERENT FROM THE "own granularity, never frame length" RULE ABOVE: a
% mode='perRunPoly' parameter's .params/.start/.lower/.upper ARE expanded to [1 x T] (one value per
% frame, the fitted polynomial evaluated at that frame's own within-run position) -- INDISTINGUISHABLE
% downstream from an ordinary 'perFrame' result, so every consumer (fitVesselPatchTimeSeries.m's own
% splitPerRun, plotVessels.m, ...) needs no special-casing at all. The RAW fitted coefficients are
% SEPARATELY exposed at .runPolyCoef.<param> (a {1 x nRunGroups} cell, each a [1 x nCoef] vector in
% .runPolyDeg.<param> order) and .runPolyDeg.<param> (the resolved degree list) -- the only place the
% actual free fit variables are visible directly. KNOWN LIMITATION carried from that expansion: .lower/
% .upper are a LINEAR expansion of a coefficient-space bound through a possibly-signed Legendre basis
% column, so the expanded number can legitimately fall outside what was actually enforced at that frame
% -- consult .runPolyCoef for the real constraint.
%
% ANGLE WRAPPING -- componentModelSpec.m's own .wrapParamPeriod (paramName -> period, e.g. 'gaussian'.
% theta -> pi, see that file) is consulted ONLY here, on the EXPANDED per-frame report: a periodic
% parameter's raw polynomial trajectory can legitimately wander outside one period as it drifts (nothing
% about the forward model cares -- cos/sin are already periodic), so .params/.start/.lower/.upper are
% wrapped into (-period/2, period/2] for reporting sanity. .runPolyCoef is NOT wrapped (it is
% coefficient space, not an angle). Every other mode is left exactly as before -- this is new machinery
% specific to the polynomial expansion, not a general per-parameter wrap applied everywhere.
%
% JOINT PATH ONLY -- there is no fitFastPath equivalent. fitVessel.m routes an all-'perFrame' fit to T
% independent small lsqnonlin calls, which is genuinely cheaper when nothing couples frames; here that
% case is close to unreachable, because an ESTIMATED 'perVessel' background (the normal way to use this
% function) already couples every frame through one shared parameter. Carrying a second engine for a
% case that essentially never fires is not worth the drift risk. The one real thing lost with
% fitFastPath is its per-frame .resnorm vector, which is given back explicitly (see .resnormPerFrame
% below) for one extra residual evaluation after convergence.
%
% LINEARITY (why summing components is exact, and why the reduction happens PER COMPONENT).
% reduceIsochromat.m is a mean over each voxel's own isochromats -- a linear operator. So
%     reduce(sum_c predIso_c) == sum_c reduce(predIso_c)
% algebraically, and this function uses the RIGHT-hand form: each component is reduced to voxel
% resolution on its own, and the reduced [nValidVox x 1] vectors are summed. Two reasons, both real:
%   - MEMORY/CACHE: a reduced contribution is nValidVox doubles, an isochromat-grid one is
%     nIso^2 times larger. Caching the reduced form (see FRAME-INVARIANT CACHING) costs almost nothing;
%     caching the grid form would not.
%   - K=1 EXACTNESS: with a single component, sum-of-one-reduction is bit-for-bit the same arithmetic
%     fitVessel.m performs, which is what makes the K=1 equivalence check meaningful rather than
%     approximate.
% The two forms can still differ in the last ulp for K>1 (floating-point addition is not associative);
% that is ordinary and not a correctness concern.
%
% FRAME-INVARIANT CACHING (opts.cacheFrameInvariant, default true) -- the optimization that makes K
% components affordable. A component whose parameters are ALL 'perVessel' (or fixed at a scalar)
% predicts the SAME image on every frame, so it need not be re-evaluated T times per residual call. Each
% component gets its own frame->group key (the unique-rows combination of its own parameters' group
% vectors), and is evaluated ONCE PER DISTINCT KEY, then reused. Cost per residual call drops from
% (nComponent * T) model evaluations to sum_c nGroups_c. For the intended usage -- a main peak free per
% frame, secondary peaks pinned per vessel, one shared background -- that is the difference between
% paying K times over and paying roughly once. This is a pure reuse of already-computed values, so it is
% BIT-IDENTICAL to recomputing (verified by A/B against opts.cacheFrameInvariant=false, not assumed);
% the option exists only so that check can be run.
%
% CONTIGUOUS GROUP LABELS ARE ENFORCED, not assumed. opts.runIdx/opts.knotIdx must be exactly
% 1:nGroups with no gaps. fitVessel.m documents a KNOWN LIMITATION here -- a grouping vector with a gap
% (e.g. [1 1 3 3]) reads past its own allocated block and crashes deep inside the residual with an
% opaque index error. Since this is a new file, that is checked up front instead
% (fitComponentFit:nonContiguousRunIdx / :nonContiguousKnotIdx), so the failure names the actual
% problem. Same "fail loud, never silently relabel" call fitVessel.m made, just earlier and legibly.
%
% NO-ARG CALL -- fitComponentFit() prints every opts field with its allowed values and returns opts
% fully populated with defaults, same convention as fitVessel.m/fitPatchVessels.m.
%
% INPUT
%   im4        : [ny x nx x 1 x T] double, ONE run-concatenated series (this function does no run
%                combination of its own -- the caller has already decided what a "frame" is).
%   validMask  : [ny x nx] logical, the voxels actually fit against. Predictions and residuals are
%                produced in the same linear order frame(validMask) gives.
%   voxSz2     : [rowMm colMm] in-plane voxel size, for the isochromat grid's own mm coordinates.
%   components : 1xC struct array, one per additive component. Each entry:
%       .model        char naming the model, OR a spec struct as returned by componentModelSpec.m.
%                     Passing a spec directly lets a caller (or a test) use a model that is not in the
%                     registry -- the engine only ever consumes the spec's fields.
%       .label        OPTIONAL char, carried through to the output untouched (provenance; e.g. which
%                     peak this component is). Default ''.
%       .mode         struct, one field per parameter: 'perVessel'|'perFrame'|'perRun'|'perKnot'|
%                     'perRunPoly([<degrees>])' (see that section above). Default 'perFrame' for any
%                     parameter left unset.
%       .fixed        struct, one field per parameter, logical. Default false. A FIXED parameter is
%                     held at its own .start value and is not a free variable -- with .start now being
%                     a literal number, that works with no seed/source machinery of any kind.
%                     mode='perRunPoly' does not support fixed=true (see that section).
%       .start        struct, one field per parameter, REQUIRED for every parameter (there are no
%                     default-heuristic start values in this engine -- see PURE NUMERICAL ENGINE).
%                     Scalar (broadcast to every group) or exactly nGroups long for that parameter's
%                     own granularity. mode='perRunPoly' additionally requires this to be SCALAR (see
%                     SEEDING above -- it seeds the degree-0 coefficient of every run).
%       .lower/.upper struct, one field per parameter. Same shape rule as .start. Default -Inf/+Inf.
%                     Ignored for a fixed parameter (not a free variable, so it has no bound).
%                     mode='perRunPoly' additionally requires these to be SCALAR (see SEEDING above).
%   opts       : struct --
%       .nIsochromatPerVoxDim  positive ODD integer -- see buildIsochromatGrid.m. Default 7.
%       .runIdx    [1 x T] frame->run group, values exactly 1:nRun. Default ones(1,T) (single run).
%       .knotIdx   [1 x T] frame->knot group, values exactly 1:nKnot. Default [] -- REQUIRED if any
%                  parameter is an ESTIMATED 'perKnot' (fitComponentFit:knotIdxRequired otherwise).
%       .cacheFrameInvariant  logical, default TRUE -- see FRAME-INVARIANT CACHING. false recomputes
%                  every component on every frame; same answer, slower. For A/B checking only.
%       .display   lsqnonlin 'Display'. Default 'off'.
%
% OUTPUT
%   out : struct --
%     .components  1xC struct array, aligned with the input. Each entry: .model/.label/.mode/.fixed
%                  (echoed), and .params/.start/.lower/.upper -- four structs, one field per parameter.
%                  Each parameter's value is reported AT ITS OWN GRANULARITY (scalar for 'perVessel',
%                  1xnRun for 'perRun', 1xnKnot for 'perKnot', 1xT for 'perFrame'), never broadcast out
%                  to frame length -- the same convention fitVessel.m's own DIAGNOSTIC OUTPUT documents.
%                  ONE DELIBERATE EXCEPTION: mode='perRunPoly' also reports 1xT (see that section's own
%                  OUTPUT note) -- indistinguishable downstream from 'perFrame', with the actual free
%                  coefficients living separately at .runPolyCoef/.runPolyDeg (present only for a
%                  component that has at least one mode='perRunPoly' parameter).
%                  A fixed parameter reports its own value for .params/.start/.lower/.upper alike.
%                  (Parameters live under .params rather than directly on the entry so that every
%                  component has the SAME outer fields regardless of model, which is what lets this be
%                  a struct array at all. A caller wanting them flat -- e.g.
%                  fitVesselPatchTimeSeries.m's own opts.model-shaped output -- reshapes at its own
%                  boundary.)
%     .runPolyCoef / .runPolyDeg  per component, one field per mode='perRunPoly' parameter (absent for
%                  every other parameter) -- see that section's own OUTPUT note.
%     .resnorm     scalar, lsqnonlin's own sum of squared residuals over the whole joint fit.
%     .resnormPerFrame  [1 x T], that same total split per frame -- what fitFastPath used to give for
%                  free and the joint path never did. Costs one extra residual evaluation. Answers
%                  "which frame did this fit fall apart on", which a per-frame multi-peak fit invites.
%     .exitFlag / .iterations / .maxIterations / .functionTolerance  lsqnonlin convergence diagnostics
%                  (iterations==maxIterations is a red flag -- stopped on the cap, not converged).
%     .nFree       number of free scalar variables actually optimized.
%
% See also componentModelSpec, fitVessel, buildIsochromatGrid, reduceIsochromat.

    if nargin==0
        printOptsHelp();
        out = defaultOpts();
        return
    end
    if nargin<5; opts = struct(); end

    %% ---- shapes, opts, validation ------------------------------------------
    im4 = double(im4(:,:,1,:));
    [ny,nx,~,T] = size(im4);
    opts = fillOptsFromDefaults(opts, defaultOpts());
    if isempty(opts.runIdx); opts.runIdx = ones(1,T); end
    runIdx = opts.runIdx(:).';
    knotIdx = opts.knotIdx(:).';
    nIso = opts.nIsochromatPerVoxDim;

    assert(isscalar(nIso) && nIso>=1 && mod(nIso,2)==1, 'fitComponentFit:nIsoMustBeOdd', ...
        'opts.nIsochromatPerVoxDim must be a positive odd integer (got %s).', mat2str(nIso));
    assert(islogical(validMask) && isequal(size(validMask),[ny nx]), 'fitComponentFit:badValidMask', ...
        'validMask must be a %dx%d logical, matching im4''s own first two dimensions.', ny, nx);
    assert(any(validMask(:)), 'fitComponentFit:emptyValidMask', 'validMask selects no voxel at all.');
    assert(numel(runIdx)==T, 'fitComponentFit:badRunIdxLength', ...
        'opts.runIdx has %d entry(ies), expected %d (one per frame).', numel(runIdx), T);
    nRunGroups = checkContiguous(runIdx, 'runIdx');
    nKnotGroups = 0;
    if ~isempty(knotIdx)
        assert(numel(knotIdx)==T, 'fitComponentFit:badKnotIdxLength', ...
            'opts.knotIdx has %d entry(ies), expected %d (one per frame).', numel(knotIdx), T);
        nKnotGroups = checkContiguous(knotIdx, 'knotIdx');
    end
    nValidVox = nnz(validMask);

    %% ---- normalize components, resolve every parameter ---------------------
    C = numel(components);
    assert(C>=1, 'fitComponentFit:noComponents', ...
        'components is empty -- at least one additive component is required.');
    % Built by assigning element 1 first rather than repmat'ing a template: every element then comes
    % from the same function and so has identical field ORDER, which struct-array assignment requires
    % (a same-fields-different-order template raises "Subscripted assignment between dissimilar
    % structures").
    comp = normalizeComponent(components(1), 1, T, runIdx, knotIdx, nRunGroups, nKnotGroups);
    for c = 2:C
        comp(c) = normalizeComponent(components(c), c, T, runIdx, knotIdx, nRunGroups, nKnotGroups);
    end

    %% ---- flat FREE-parameter slot list, ordered by granularity class -------
    % One slot per FREE (component, parameter) pair. Ordered perVessel, then perFrame, then perRun,
    % then perKnot, then perRunPoly LAST -- and within each class by component order, then that model's
    % own paramNames order. The ordering is deliberate rather than incidental: with a single component
    % (and no perRunPoly parameter) it reproduces fitVessel.m's own free-vector layout exactly, which is
    % what allows the K=1 equivalence check to expect bit-identical output rather than merely close
    % output. perRunPoly is appended in its own pass (below) rather than folded into modeOrder because
    % its canonical mode STRING varies per parameter (the embedded degree list) -- P.isPoly, not a
    % strcmp against one literal, is what selects it.
    modeOrder = {'perVessel','perFrame','perRun','perKnot'};
    slots = struct('comp',{}, 'param',{}, 'mode',{}, 'groupIdx',{}, 'nGroups',{}, 'offset',{});
    for mo = 1:numel(modeOrder)
        for c = 1:C
            for p = 1:numel(comp(c).paramNames)
                nm = comp(c).paramNames{p};
                P = comp(c).param.(nm);
                if P.fixed || P.isPoly || ~strcmp(P.mode, modeOrder{mo}); continue; end
                slots(end+1) = mkSlot(c, nm, P); %#ok<AGROW>
            end
        end
    end
    for c = 1:C
        for p = 1:numel(comp(c).paramNames)
            nm = comp(c).paramNames{p};
            P = comp(c).param.(nm);
            if P.fixed || ~P.isPoly; continue; end
            slots(end+1) = mkSlot(c, nm, P); %#ok<AGROW>
        end
    end
    nSlot = numel(slots);
    assert(nSlot>0, 'fitComponentFit:nothingToFit', ...
        'every parameter of every component is fixed -- there is nothing to optimize.');

    nFree = 0;
    for s = 1:nSlot
        slots(s).offset = nFree;          % 0-based start of this slot's own column block
        nFree = nFree + slots(s).nGroups;
        % Back-reference so the residual can go parameter -> slot directly, instead of searching the
        % slot list once per parameter per frame per lsqnonlin call.
        comp(slots(s).comp).param.(slots(s).param).slotIdx = s;
    end
    assert(nValidVox*T >= nFree, 'fitComponentFit:underdetermined', ...
        ['%d free variable(s) against only %d residual(s) (%d valid voxel(s) x %d frame(s)) -- the ' ...
         'fit is underdetermined. Coarsen a granularity, fix a parameter, or widen validMask.'], ...
        nFree, nValidVox*T, nValidVox, T);

    %% ---- start point and bounds -------------------------------------------
    x0v = zeros(1,nFree); lbv = zeros(1,nFree); ubv = zeros(1,nFree);
    for s = 1:nSlot
        P = comp(slots(s).comp).param.(slots(s).param);
        idx = slots(s).offset + (1:slots(s).nGroups);
        x0v(idx) = P.start;
        lbv(idx) = P.lower;
        ubv(idx) = P.upper;
    end

    %% ---- data, isochromat grid, per-component evaluators -------------------
    ctrCol = (nx+1)/2; ctrRow = (ny+1)/2;
    [XIso, YIso] = buildIsochromatGrid(nIso, nx, ny, ctrCol, ctrRow, voxSz2);
    zvAll = zeros(nValidVox, T);
    for t = 1:T
        frame = im4(:,:,1,t);
        nanM = isnan(frame);
        if any(nanM(:)); frame(nanM) = mean(frame(~nanM),'omitnan'); end
        zvAll(:,t) = frame(validMask);
    end

    S = struct();
    S.comp = comp; S.C = C; S.slots = slots; S.nSlot = nSlot; S.T = T;
    S.XIso = XIso; S.YIso = YIso; S.nIso = nIso; S.ny = ny; S.nx = nx;
    S.validMask = validMask; S.zvAll = zvAll; S.nValidVox = nValidVox;
    S.cache = opts.cacheFrameInvariant;

    %% ---- solve --------------------------------------------------------------
    % MaxFunctionEvaluations scales with the SLOT count (free parameters), not the column count -- the
    % same quantity fitVessel.m scales by, so a K=1 fit gets the identical budget.
    optsLs = optimoptions('lsqnonlin', 'Display',opts.display, ...
        'MaxFunctionEvaluations', 3000*max(1,nSlot), ...
        'JacobPattern', buildSlotJacobPattern(slots, comp, nValidVox, T, nFree));
    residFun = @(xx) evalResidualComponents(xx, S);
    % SIX AND SEVEN MATTER. lambda carries the BOUND multipliers and jacobian the Jacobian at the
    % solution; both were previously discarded, and both are what makes parameter identifiability
    % answerable from the solve itself rather than from a size heuristic. See IDENTIFIABILITY below
    % and fitIdentifiability.m, which is the consumer.
    [xFit, resnorm, residFinal, exitFlag, lsqOutput, lsqLambda, lsqJacobian] = ...
        lsqnonlin(residFun, x0v, lbv, ubv, optsLs);

    %% ---- unpack -------------------------------------------------------------
    outComp = repmat(struct('model','', 'label','', 'mode',struct(), 'fixed',struct(), ...
                            'params',struct(), 'start',struct(), 'lower',struct(), 'upper',struct(), ...
                            'runPolyCoef',struct(), 'runPolyDeg',struct()), 1, C);
    for c = 1:C
        outComp(c).model = comp(c).spec.name;
        outComp(c).label = comp(c).label;
        outComp(c).mode  = comp(c).mode;
        outComp(c).fixed = comp(c).fixed;
        for p = 1:numel(comp(c).paramNames)
            nm = comp(c).paramNames{p};
            P = comp(c).param.(nm);
            if P.fixed
                % Reported at its own granularity, NOT broadcast to frame length -- the convention
                % fitVessel.m documents and (as of the stage-0 fix) honours. start/lower/upper all
                % report the same value: a fixed parameter has no real bound.
                outComp(c).params.(nm) = P.startRaw;
                outComp(c).start.(nm)  = P.startRaw;
                outComp(c).lower.(nm)  = P.startRaw;
                outComp(c).upper.(nm)  = P.startRaw;
            elseif P.isPoly
                % See file header mode='perRunPoly' OUTPUT note: expanded to [1 x T] (indistinguishable
                % downstream from 'perFrame'), raw coefficients kept separately at .runPolyCoef/
                % .runPolyDeg, and -- for a genuinely periodic parameter (componentModelSpec.m's own
                % .wrapParamPeriod) -- wrapped into its principal range.
                sIdx = findSlot(slots, c, nm);
                off = slots(sIdx).offset; nCoefK = P.nCoefPerRun; nRunG = numel(P.basisByRun);
                coefFit = reshape(xFit(off+(1:P.nGroups)), nCoefK, nRunG).';
                coefSt  = reshape(x0v(off+(1:P.nGroups)),  nCoefK, nRunG).';
                coefLo  = reshape(lbv(off+(1:P.nGroups)),  nCoefK, nRunG).';
                coefUp  = reshape(ubv(off+(1:P.nGroups)),  nCoefK, nRunG).';
                vFit = zeros(1,T); vSt = vFit; vLo = vFit; vUp = vFit;
                for g = 1:nRunG
                    idxG = P.runIdxLocal == g;
                    Bg = P.basisByRun{g};
                    vFit(idxG) = (Bg * coefFit(g,:).').';
                    vSt(idxG)  = (Bg * coefSt(g,:).').';
                    vLo(idxG)  = (Bg * coefLo(g,:).').';
                    vUp(idxG)  = (Bg * coefUp(g,:).').';
                end
                period = wrapPeriodFor(comp(c).spec, nm);
                if ~isempty(period)
                    vFit = wrapAngleLocal(vFit, period);
                    vSt  = wrapAngleLocal(vSt,  period);
                    vLo  = wrapAngleLocal(vLo,  period);
                    vUp  = wrapAngleLocal(vUp,  period);
                end
                outComp(c).params.(nm) = vFit;
                outComp(c).start.(nm)  = vSt;
                outComp(c).lower.(nm)  = vLo;
                outComp(c).upper.(nm)  = vUp;
                outComp(c).runPolyCoef.(nm) = arrayfun(@(g) coefFit(g,:), 1:nRunG, 'UniformOutput', false);
                outComp(c).runPolyDeg.(nm)  = P.polyDegs;
            else
                sIdx = findSlot(slots, c, nm);
                idx = slots(sIdx).offset + (1:slots(sIdx).nGroups);
                outComp(c).params.(nm) = xFit(idx);
                outComp(c).start.(nm)  = x0v(idx);
                outComp(c).lower.(nm)  = lbv(idx);
                outComp(c).upper.(nm)  = ubv(idx);
            end
        end
        outComp(c) = canonicalizeComponent(outComp(c), comp(c));
    end

    out = struct();
    out.components        = outComp;
    out.resnorm           = resnorm;
    out.resnormPerFrame   = sum(reshape(residFinal, nValidVox, T).^2, 1);
    out.exitFlag          = exitFlag;
    out.iterations        = lsqOutput.iterations;
    out.maxIterations     = optsLs.MaxIterations;
    out.functionTolerance = optsLs.FunctionTolerance;
    out.nFree             = nFree;

    % IDENTIFIABILITY -- RAW solver output only, nothing derived from it. Same rule
    % fitVesselPatchTimeSeries.m applies to replicates: the producer produces, the consumer derives.
    % fitIdentifiability.m turns these into a spectrum / null directions / active-bound flags.
    %   .jacobian   [m x nFree] SPARSE, at the solution. Finite-difference under JacobPattern (no
    %               analytic jacobian is supplied), so it carries FD step error and the imposed
    %               sparsity -- fine for conditioning, not exact.
    %   .lambda     lsqnonlin's bound multipliers (.lower/.upper, nFree each). Carried because it is
    %               raw solver output and costs nothing, but NOT the way to find the active set:
    %               these are not zero in the interior, so testing them against zero asks "does this
    %               parameter have a bound?" rather than "does it rest on one?". Measured: that test
    %               flagged 16 of 19 parameters on every real vessel, sparing only the unbounded
    %               theta columns. Use |p-bound| instead -- .freeLower/.freeUpper below exist for it.
    %               The active set still has to be found somehow: for a parameter resting on a bound
    %               the jacobian says nothing useful about its uncertainty, so a J'J error bar there
    %               is meaningless rather than merely imprecise.
    %   .freeLabel  1 x nFree cellstr naming each free column '<label>.<param>[group g]', so a null
    %               DIRECTION can be reported in terms a reader recognises instead of column indices.
    %   .freeComp / .freeParam / .freeGroup  the same mapping in machine-readable form.
    %   .freeStart/.freeLower/.freeUpper  the solver's own x0/lb/ub in FREE-VECTOR order. The bounds
    %               are what let a consumer test |p-bound| directly instead of trusting the
    %               multipliers alone -- a multiplier can be exactly zero at a bound the solver only
    %               touched, and a missed active bound silently yields a confident error bar for a
    %               number the data never set.
    out.jacobian  = lsqJacobian;
    out.lambda    = lsqLambda;
    out.freeStart = x0v(:).';
    out.freeLower = lbv(:).';
    out.freeUpper = ubv(:).';
    out.freeFit   = xFit(:).';
    [out.freeLabel, out.freeComp, out.freeParam, out.freeGroup] = buildFreeLabels(slots, outComp);
end

% ---------------------------------------------------------------------------
function [lab, cIdx, pName, gIdx] = buildFreeLabels(slots, outComp)
    % One entry per FREE column, in free-vector order. Built from slots (the packing authority) so it
    % cannot drift from the actual column layout.
    nF = 0;
    for s = 1:numel(slots); nF = nF + slots(s).nGroups; end
    lab = cell(1, nF); cIdx = zeros(1, nF); pName = cell(1, nF); gIdx = zeros(1, nF);
    for s = 1:numel(slots)
        c = slots(s).comp; nm = slots(s).param;
        for g = 1:slots(s).nGroups
            k = slots(s).offset + g;
            if slots(s).nGroups == 1
                lab{k} = sprintf('%s.%s', outComp(c).label, nm);
            else
                lab{k} = sprintf('%s.%s[%s %d]', outComp(c).label, nm, slots(s).mode, g);
            end
            cIdx(k) = c; pName{k} = nm; gIdx(k) = g;
        end
    end
end

% ---------------------------------------------------------------------------
function nGroups = checkContiguous(idxVec, name)
    % Group labels must be exactly 1:nGroups. See file header CONTIGUOUS GROUP LABELS.
    u = unique(idxVec);
    nGroups = numel(u);
    assert(isequal(u(:).', 1:nGroups), ['fitComponentFit:nonContiguous' upper(name(1)) name(2:end)], ...
        ['opts.%s must label its groups exactly 1:%d with no gaps (got unique values %s). A gap would ' ...
         'index past a parameter''s own allocated column block; relabel the groups contiguously rather ' ...
         'than leaving an unused label.'], name, nGroups, mat2str(u(:).'));
end

% ---------------------------------------------------------------------------
function cOut = normalizeComponent(cIn, cIdx, T, runIdx, knotIdx, nRunGroups, nKnotGroups)
    % Resolve one component: its spec, every parameter's mode/fixed/start/lower/upper, and the
    % frame->group vector that IS its granularity (see file header).
    assert(isfield(cIn,'model') && ~isempty(cIn.model), 'fitComponentFit:componentNoModel', ...
        'components(%d) has no .model (a model name, or a componentModelSpec.m spec struct).', cIdx);
    if isstruct(cIn.model); spec = cIn.model; else; spec = componentModelSpec(cIn.model); end
    assert(all(isfield(spec, {'name','paramNames','modelFun','toModelParams'})), ...
        'fitComponentFit:badModelSpec', ...
        'components(%d).model resolved to a struct missing one of .name/.paramNames/.modelFun/.toModelParams.', cIdx);

    cOut = struct();
    cOut.spec = spec;
    cOut.paramNames = spec.paramNames;
    cOut.label = '';
    if isfield(cIn,'label') && ~isempty(cIn.label); cOut.label = char(string(cIn.label)); end
    cOut.mode = struct(); cOut.fixed = struct(); cOut.param = struct();

    grpAll = cell(1, numel(spec.paramNames));
    for p = 1:numel(spec.paramNames)
        nm = spec.paramNames{p};
        P = struct();
        P.mode  = getSub(cIn, 'mode',  nm, 'perFrame');
        P.fixed = logical(getSub(cIn, 'fixed', nm, false));
        P.isPoly = false;

        [isPolyMode, polyDegs, canonicalMode] = parsePerRunPolyMode(P.mode);
        if isPolyMode
            P.mode = canonicalMode;
        else
            assert(ismember(P.mode, {'perVessel','perFrame','perRun','perKnot'}), ...
                'fitComponentFit:badMode', ...
                ['components(%d) (%s) parameter ''%s'': unknown mode ''%s'' (perVessel/perFrame/perRun/' ...
                 'perKnot/perRunPoly([<degrees>])).'], cIdx, spec.name, nm, char(string(P.mode)));
        end

        if isPolyMode
            % See file header mode='perRunPoly'. A design MATRIX, not a group selector, so it takes its
            % own path rather than the frame->group switch below.
            assert(~P.fixed, 'fitComponentFit:perRunPolyFixedUnsupported', ...
                ['components(%d) (%s) parameter ''%s'': mode=''perRunPoly'' does not support fixed=true ' ...
                 '-- pin a run''s own within-run trend by fitting/precomputing it separately per run ' ...
                 'instead (mode=''perRun'', fixed=true, literal per-run .start), or leave this parameter ' ...
                 'estimated.'], cIdx, spec.name, nm);
            nCoef = numel(polyDegs);
            runLensLocal = arrayfun(@(g) nnz(runIdx==g), 1:nRunGroups);
            assert(all(nCoef < runLensLocal), 'fitComponentFit:runPolyDegSaturated', ...
                ['components(%d) (%s) parameter ''%s'': %d within-run polynomial term(s) (degrees %s) ' ...
                 'requested but run(s) have as few as %d frame(s) -- a saturated/degenerate basis.'], ...
                cIdx, spec.name, nm, nCoef, mat2str(polyDegs), min(runLensLocal));
            localFrameIdx = zeros(1,T);
            basisByRun = cell(1,nRunGroups);
            for g = 1:nRunGroups
                idxG = find(runIdx==g);
                localFrameIdx(idxG) = 1:numel(idxG);
                basisByRun{g} = legendreBasisLocal(runLensLocal(g), polyDegs);
            end
            P.isPoly = true;
            P.polyDegs = polyDegs;
            P.nCoefPerRun = nCoef;
            P.basisByRun = basisByRun;
            P.runIdxLocal = runIdx;
            P.localFrameIdx = localFrameIdx;
            % CACHING-KEY ONLY (see FRAME-INVARIANT CACHING) -- NOT used for direct free-vector
            % indexing (evalResidualComponents branches on P.isPoly before ever consulting groupIdx for
            % a poly parameter). Every frame generally differs under a polynomial trajectory, so this
            % is exactly 'perFrame''s own contribution: no false caching, and none missed either.
            P.groupIdx = 1:T;
            P.nGroups = nCoef * nRunGroups;

            startScalar = getSub(cIn, 'start', nm, []);
            assert(~isempty(startScalar) && isnumeric(startScalar) && isscalar(startScalar), ...
                'fitComponentFit:perRunPolyStartMustBeScalar', ...
                ['components(%d) (%s) parameter ''%s'': mode=''perRunPoly'' requires a SCALAR .start -- ' ...
                 'it seeds the degree-0 coefficient of every run (every other coefficient starts at 0); ' ...
                 'got %d value(s).'], cIdx, spec.name, nm, numel(startScalar));
            P.startRaw = startScalar;
            lowerScalar = getSub(cIn, 'lower', nm, -Inf);
            upperScalar = getSub(cIn, 'upper', nm,  Inf);
            assert(isscalar(lowerScalar) && isscalar(upperScalar), ...
                'fitComponentFit:perRunPolyBoundsMustBeScalar', ...
                ['components(%d) (%s) parameter ''%s'': mode=''perRunPoly'' requires SCALAR .lower/' ...
                 '.upper -- applied to the degree-0 coefficient of every run; every other coefficient ' ...
                 'gets a default symmetric +-(upper-lower)/2 bound.'], cIdx, spec.name, nm);
            halfWidth = (upperScalar - lowerScalar) / 2;
            startVec = zeros(1, P.nGroups); lowerVec = startVec; upperVec = startVec;
            for g = 1:nRunGroups
                for c = 1:nCoef
                    k = (g-1)*nCoef + c;
                    if polyDegs(c) == 0
                        startVec(k) = startScalar; lowerVec(k) = lowerScalar; upperVec(k) = upperScalar;
                    else
                        startVec(k) = 0; lowerVec(k) = -halfWidth; upperVec(k) = halfWidth;
                    end
                end
            end
            P.start = startVec; P.lower = lowerVec; P.upper = upperVec;
            tolB = 1e-9 * max(1, abs(P.start));
            assert(all(P.start >= P.lower - tolB) && all(P.start <= P.upper + tolB), ...
                'fitComponentFit:startOutsideBounds', ...
                ['components(%d) (%s) parameter ''%s'': .start falls outside [.lower .upper] once ' ...
                 'expanded per run/coefficient (start %s, lower %s, upper %s).'], ...
                cIdx, spec.name, nm, mat2str(P.start,4), mat2str(P.lower,4), mat2str(P.upper,4));
        else
            % The frame->group vector, and with it the group count -- this IS the granularity.
            switch P.mode
                case 'perVessel'; P.groupIdx = ones(1,T);  P.nGroups = 1;
                case 'perFrame';  P.groupIdx = 1:T;        P.nGroups = T;
                case 'perRun';    P.groupIdx = runIdx;     P.nGroups = nRunGroups;
                case 'perKnot'
                    assert(~isempty(knotIdx), 'fitComponentFit:knotIdxRequired', ...
                        ['components(%d) (%s) parameter ''%s'' has mode=''perKnot'' but opts.knotIdx is ' ...
                         'empty -- a [1 x T] frame->knot-group vector is required (this engine has no ' ...
                         'knot concept of its own).'], cIdx, spec.name, nm);
                    P.groupIdx = knotIdx; P.nGroups = nKnotGroups;
            end

            startRaw = getSub(cIn, 'start', nm, []);
            assert(~isempty(startRaw) && isnumeric(startRaw), 'fitComponentFit:startRequired', ...
                ['components(%d) (%s) parameter ''%s'': .start is REQUIRED and must be numeric -- this ' ...
                 'engine has no default start-value heuristic (see file header PURE NUMERICAL ENGINE); ' ...
                 'the caller resolves seeding.'], cIdx, spec.name, nm);
            P.startRaw = startRaw(:).';
            P.start = fitToGroups(P.startRaw, P.nGroups, cIdx, spec.name, nm, 'start');
            if P.fixed
                % Expanded onto every frame, which is what the residual indexes. The UN-expanded
                % .startRaw is kept for reporting -- the same two-shapes split fitVessel.m uses
                % (.fixedVal/.outVal).
                P.fixedVal = P.start(P.groupIdx);
            else
                P.lower = fitToGroups(getSub(cIn,'lower',nm,-Inf), P.nGroups, cIdx, spec.name, nm, 'lower');
                P.upper = fitToGroups(getSub(cIn,'upper',nm, Inf), P.nGroups, cIdx, spec.name, nm, 'upper');
                % Tolerance derived from START, never from the BOUND: a bound of +/-Inf is the normal
                % case here (theta is always fit completely unbounded in this pipeline -- see
                % fitVessel.m's own THETA CANONICALIZATION), and eps(Inf) is NaN, which would silently
                % poison the comparison and reject every real call. Inf +/- a finite tolerance stays
                % Inf, so this form is safe.
                tolB = 1e-9 * max(1, abs(P.start));
                assert(all(P.start >= P.lower - tolB) && all(P.start <= P.upper + tolB), ...
                    'fitComponentFit:startOutsideBounds', ...
                    ['components(%d) (%s) parameter ''%s'': .start falls outside [.lower .upper] ' ...
                     '(start %s, lower %s, upper %s) -- lsqnonlin would silently clip it.'], ...
                    cIdx, spec.name, nm, mat2str(P.start,4), mat2str(P.lower,4), mat2str(P.upper,4));
            end
        end

        cOut.mode.(nm) = P.mode; cOut.fixed.(nm) = P.fixed;
        cOut.param.(nm) = P;
        grpAll{p} = P.groupIdx;
    end

    % This component's OWN frame->group key: the unique combination of its parameters' group vectors.
    % All-'perVessel' collapses to a single key (evaluate once per residual call); any 'perFrame'
    % parameter forces T. 'stable' so key 1 is frame 1's group -- deterministic and readable.
    G = cat(1, grpAll{:});
    [~, ~, keyIdx] = unique(G.', 'rows', 'stable');
    cOut.keyIdx = keyIdx(:).';
    cOut.nKey = max(cOut.keyIdx);
    cOut.repFrame = zeros(1, cOut.nKey);
    for u = 1:cOut.nKey; cOut.repFrame(u) = find(cOut.keyIdx==u, 1); end
end

% ---------------------------------------------------------------------------
function v = getSub(cIn, group, name, dflt)
    v = dflt;
    if isfield(cIn, group) && isstruct(cIn.(group)) && isfield(cIn.(group), name) ...
            && ~isempty(cIn.(group).(name))
        v = cIn.(group).(name);
    end
end

% ---------------------------------------------------------------------------
function v = fitToGroups(v, nGroups, cIdx, modelName, paramName, what)
    % Scalar broadcasts to every group; otherwise the length must match the granularity exactly.
    v = v(:).';
    if isscalar(v)
        v = repmat(v, 1, nGroups);
    else
        assert(numel(v)==nGroups, 'fitComponentFit:badValueLength', ...
            ['components(%d) (%s) parameter ''%s'': .%s has %d value(s), expected 1 (broadcast) or %d ' ...
             '(one per group at this parameter''s own granularity).'], ...
            cIdx, modelName, paramName, what, numel(v), nGroups);
    end
end

% ---------------------------------------------------------------------------
function sIdx = findSlot(slots, c, name)
    sIdx = find([slots.comp]==c & strcmp({slots.param}, name), 1);
end

% ---------------------------------------------------------------------------
function res = evalResidualComponents(xx, S)
    % pred(t) = sum_c reduce(modelFun_c(...)) -- reduced PER COMPONENT then summed (see file header
    % LINEARITY), with each component evaluated once per distinct frame-group key rather than once per
    % frame (see FRAME-INVARIANT CACHING).
    predAll = zeros(S.nValidVox, S.T);
    for c = 1:S.C
        cc = S.comp(c);
        if S.cache
            keys = 1:cc.nKey; frames = cc.repFrame;   % one evaluation per distinct key
        else
            keys = 1:S.T;     frames = 1:S.T;         % A/B mode: recompute every frame
        end
        red = zeros(S.nValidVox, numel(keys));
        for k = 1:numel(keys)
            t = frames(k);
            vals = cell(1, numel(cc.paramNames));
            for p = 1:numel(cc.paramNames)
                nm = cc.paramNames{p};
                P = cc.param.(nm);
                if P.fixed
                    vals{p} = P.fixedVal(t);
                elseif P.isPoly
                    % Design-MATRIX value, not a group lookup -- see file header mode='perRunPoly'.
                    % localFr, not k -- k is the ENCLOSING loop's own key index (see red(:,k) below);
                    % shadowing it here silently corrupted that indexing (caught by checkcode).
                    off = S.slots(P.slotIdx).offset;
                    g = P.runIdxLocal(t); localFr = P.localFrameIdx(t); nCoefK = P.nCoefPerRun;
                    coefIdx = off + (g-1)*nCoefK + (1:nCoefK);
                    vals{p} = P.basisByRun{g}(localFr,:) * xx(coefIdx).';
                else
                    sIdx = P.slotIdx;
                    vals{p} = xx(S.slots(sIdx).offset + S.slots(sIdx).groupIdx(t));
                end
            end
            modelVals = cc.spec.toModelParams(vals);
            predIso = cc.spec.modelFun(S.XIso, S.YIso, modelVals{:});
            red(:,k) = reduceIsochromat(predIso, S.nIso, S.ny, S.nx, S.validMask);
        end
        if S.cache
            predAll = predAll + red(:, cc.keyIdx);   % scatter each key's own value onto its frames
        else
            predAll = predAll + red;
        end
    end
    res = predAll(:) - S.zvAll(:);
end

% ---------------------------------------------------------------------------
function Jpat = buildSlotJacobPattern(slots, comp, nValidVox, T, nFree)
    % Sparsity pattern for the stacked residual: T frame-blocks of nValidVox rows each, in the order
    % evalResidualComponents stacks them (predAll(:) is column-major, i.e. frame-major). One uniform
    % rule covers every GROUP-SELECTOR granularity: slot s's column g touches exactly the frames where
    % its own groupIdx==g. A 'perVessel' slot has one column touching all frames (dense); a 'perFrame'
    % slot has T columns each touching one frame (block-diagonal); 'perRun'/'perKnot' sit in between.
    % fitVessel.m needs four separate arms for this; representing granularity as a frame->group vector
    % collapses them into the loop below. mode='perRunPoly' (P.isPoly) is NOT a group selector -- every
    % one of its nCoefPerRun columns for run g touches ALL of run g's own frames (dense within the run,
    % zero outside it, the SAME "mutually exclusive across runs" pattern 'perRun' already has, just
    % nCoefPerRun columns per run instead of 1) -- handled in its own arm, reading the poly fields off
    % comp(...).param(...) directly (slots itself carries no poly-specific fields).
    %
    % WHY IT MATTERS (unchanged from fitVessel.m's own reasoning): without a pattern, lsqnonlin's
    % finite-difference Jacobian treats every column as dense and perturbs them one at a time, and each
    % perturbation re-evaluates ALL T frames -- quadratic in T as soon as anything is 'perFrame'. With
    % the pattern, column-coloring groups every mutually-non-overlapping column into ONE perturbation,
    % so the call count is O(number of slots [or poly coefficients]), independent of T.
    m = nValidVox*T;
    frameRows = reshape(1:m, nValidVox, T);
    rowParts = {}; colParts = {};
    for s = 1:numel(slots)
        P = comp(slots(s).comp).param.(slots(s).param);
        if P.isPoly
            nCoefK = P.nCoefPerRun;
            for g = 1:numel(P.basisByRun)
                rows = frameRows(:, P.runIdxLocal==g);
                for cIdxLocal = 1:nCoefK
                    rowParts{end+1} = rows(:); %#ok<AGROW>
                    colParts{end+1} = repmat(slots(s).offset+(g-1)*nCoefK+cIdxLocal, numel(rows), 1); %#ok<AGROW>
                end
            end
        elseif slots(s).nGroups == 1
            rowParts{end+1} = (1:m).'; %#ok<AGROW>
            colParts{end+1} = repmat(slots(s).offset+1, m, 1); %#ok<AGROW>
        else
            % Columns are disjoint by construction (each frame belongs to exactly one group), so the
            % whole slot's pattern is one scatter: every frame's rows against that frame's own column.
            gi = slots(s).groupIdx;
            rowParts{end+1} = frameRows(:); %#ok<AGROW>
            colParts{end+1} = slots(s).offset + repelem(gi(:), nValidVox); %#ok<AGROW>
        end
    end
    Jpat = sparse(cat(1,rowParts{:}), cat(1,colParts{:}), 1, m, nFree);
end

% ---------------------------------------------------------------------------
function s = mkSlot(c, nm, P)
    % One entry of the flat FREE-parameter slot list (see that section's own file header note). Direct
    % field assignment, not struct(...), specifically so a poly parameter's P.basisByRun (a CELL) lands
    % in s.basisByRun literally rather than triggering struct()'s cell-expands-to-struct-array pitfall --
    % moot now that slots itself carries no poly-specific fields (see buildSlotJacobPattern), but this
    % stays the safe idiom regardless of what P carries.
    s.comp = c; s.param = nm; s.mode = P.mode;
    s.groupIdx = P.groupIdx; s.nGroups = P.nGroups; s.offset = 0;
end

% ---------------------------------------------------------------------------
function outC = canonicalizeComponent(outC, cc)
    % Per-component theta canonicalization, gated EXACTLY as fitVessel.m's own canonicalizeThetaIfFree:
    % only when both shape parameters AND theta are ESTIMATED with the IDENTICAL mode. That is the one
    % case where the shape swap is a genuine no-op, because all three were equally free to have
    % converged either way. A caller-pinned (fixed) shape parameter is left alone -- swapping it would
    % change what was actually fit, not merely relabel it.
    spec = cc.spec;
    if isempty(spec.shapeParamNames) || ~isfield(outC.params,'theta'); return; end
    s1 = spec.shapeParamNames{1}; s2 = spec.shapeParamNames{2};
    if outC.fixed.(s1) || outC.fixed.(s2) || outC.fixed.theta; return; end
    if ~strcmp(outC.mode.(s1), outC.mode.(s2)) || ~strcmp(outC.mode.(s1), outC.mode.theta); return; end
    swapIdx = spec.shapeSwapIdxFun(outC.params.(s1), outC.params.(s2));
    [outC.params.(s1), outC.params.(s2), outC.params.theta] = ...
        spec.canonicalizeShapeFun(outC.params.(s1), outC.params.(s2), outC.params.theta);
    for bf = {'start','lower','upper'}
        f = bf{1};
        [outC.(f).(s1), outC.(f).(s2)] = spec.canonicalizeShapeBoundsFun(outC.(f).(s1), outC.(f).(s2), swapIdx);
    end
end

% ---------------------------------------------------------------------------
function opts = defaultOpts()
    opts.nIsochromatPerVoxDim = 7;
    opts.runIdx  = [];   % empty placeholder -> ones(1,T) at call time (T isn't known here)
    opts.knotIdx = [];
    opts.cacheFrameInvariant = true;
    opts.display = 'off';
end

% ---------------------------------------------------------------------------
function printOptsHelp()
    fprintf('fitComponentFit opts -- allowed value(s) per field (default in parentheses):\n');
    fprintf('  opts.nIsochromatPerVoxDim : positive ODD integer                                    (7)\n');
    fprintf('  opts.runIdx               : [1 x T] frame->run group, exactly 1:nRun                ([] = single run)\n');
    fprintf('  opts.knotIdx              : [1 x T] frame->knot group, exactly 1:nKnot; REQUIRED for an estimated mode=''perKnot'' ([])\n');
    fprintf('  opts.cacheFrameInvariant  : true | false -- evaluate each component once per distinct frame-group (true)\n');
    fprintf('  opts.display              : lsqnonlin ''Display''                                     (''off'')\n');
    fprintf('components(c) fields -- .model (name or componentModelSpec struct), .label, and the\n');
    fprintf('  per-parameter structs .mode (''perVessel''|''perFrame''|''perRun''|''perKnot''|\n');
    fprintf('  ''perRunPoly([<degrees>])'', default ''perFrame'' -- see file header mode=''perRunPoly''),\n');
    fprintf('  .fixed (logical, default false; ''perRunPoly'' only supports false), .start (REQUIRED\n');
    fprintf('  numeric, SCALAR for ''perRunPoly''), .lower/.upper (numeric, default -+Inf, SCALAR for\n');
    fprintf('  ''perRunPoly'').\n');
end

% ---------------------------------------------------------------------------
function [isPoly, degs, canonicalMode] = parsePerRunPolyMode(modeStr)
    % Recognize and parse mode=''perRunPoly(<degrees>)'' -- see file header for the full syntax and the
    % four accepted spellings. Not a poly mode at all -> isPoly=false, degs/canonicalMode unset (the
    % caller falls through to the ordinary perVessel/perFrame/perRun/perKnot validation).
    isPoly = false; degs = []; canonicalMode = '';
    s = strtrim(char(string(modeStr)));
    prefix = 'perRunPoly';
    if ~strncmp(s, prefix, numel(prefix)); return; end
    isPoly = true;
    rest = strtrim(s(numel(prefix)+1:end));
    % Strip alternating layers of matching ()/[] wrapping -- 'perRunPoly([1 2])' peels twice (parens,
    % then brackets) down to '1 2'; 'perRunPoly(1,3)' peels once down to '1,3'.
    while numel(rest) >= 2 && ((rest(1)=='(' && rest(end)==')') || (rest(1)=='[' && rest(end)==']'))
        rest = strtrim(rest(2:end-1));
    end
    assert(~isempty(rest), 'fitComponentFit:perRunPolyNoTerms', ...
        ['mode ''%s'' -- mode=''perRunPoly'' requires an explicit non-negative-integer degree list, ' ...
         'e.g. ''perRunPoly([0 1 2])'', ''perRunPoly(1,3)'', ''perRunPoly[1 2]'', ''perRunPoly(0:4)''.'], s);
    tokens = regexp(rest, '[,\s]+', 'split');
    tokens = tokens(~cellfun(@isempty, tokens));
    vals = [];
    for i = 1:numel(tokens)
        v = str2num(tokens{i}); %#ok<ST2NM> -- accepts '2', '0:4', '0:2:8' (SAME idiom this codebase's
                                 % own parseSeedPath uses for its own bracketed-index selector, see
                                 % fitVesselPatchTimeSeries.m)
        assert(~isempty(v) && isnumeric(v) && all(isfinite(v)) && all(v==round(v)) && all(v>=0), ...
            'fitComponentFit:badPerRunPolyTerm', ...
            'mode ''%s'' -- term ''%s'' is not a non-negative integer (or range like ''0:4'').', s, tokens{i});
        vals = [vals v(:).']; %#ok<AGROW>
    end
    degs = unique(vals);
    canonicalMode = sprintf('perRunPoly([%s])', ...
        strjoin(arrayfun(@(d) sprintf('%d',d), degs, 'UniformOutput',false), ' '));
end

% ---------------------------------------------------------------------------
function B = legendreBasisLocal(nIdx, degs)
    % Per-run Legendre basis columns for mode='perRunPoly', for the requested DEGREES (degree 0 =
    % constant, degree>=1 DEMEANED over the run -- orthogonal to the constant, so a degree>=1
    % coefficient of 0 means exactly "no deviation from the run's own degree-0 level"). Index position
    % j=0..nIdx-1 -> x = 2j/(nIdx-1) - 1; Bonnet recurrence. Column order = sorted degs (degs is already
    % sorted/unique by the time this is called -- see parsePerRunPolyMode).
    %
    % SAME formula as fitIRF.m's own private polyBaseline/resolvePolyTerms (that file's within-run
    % baseline-drift regressor) -- duplicated rather than shared/promoted, which is out of scope for
    % this change; this codebase's established small-helper-duplication convention (see
    % componentModelSpec.m's own canonicalizeRadiusAspectShape for precedent). nIdx>=2 is guaranteed by
    % the runPolyDegSaturated assert at the one call site (normalizeComponent), before this runs -- a
    % single-frame run would otherwise divide by zero building x.
    degs = degs(:).';
    Pmax = max(degs);
    x = 2*(0:nIdx-1).'/(nIdx-1) - 1;
    L = zeros(nIdx, Pmax+1); L(:,1) = 1;
    if Pmax >= 1; L(:,2) = x; end
    for p = 1:Pmax-1; L(:,p+2) = ((2*p+1)*x.*L(:,p+1) - p*L(:,p))/(p+1); end
    L(:,2:end) = L(:,2:end) - mean(L(:,2:end), 1);
    B = L(:, degs+1);
end

% ---------------------------------------------------------------------------
function period = wrapPeriodFor(spec, nm)
    % componentModelSpec.m's own .wrapParamPeriod (paramName -> period), consulted ONLY when expanding a
    % mode='perRunPoly' parameter's fitted trajectory (see file header ANGLE WRAPPING). [] = not
    % periodic, no wrap.
    period = [];
    if isfield(spec,'wrapParamPeriod') && isstruct(spec.wrapParamPeriod) && isfield(spec.wrapParamPeriod, nm)
        period = spec.wrapParamPeriod.(nm);
    end
end

% ---------------------------------------------------------------------------
function v = wrapAngleLocal(v, period)
    % Wrap into (-period/2, period/2] -- the SAME formula this codebase already uses for theta's own
    % period-pi ellipse symmetry (componentModelSpec.m's canonicalizeRadiusAspectShape), generalized over
    % an arbitrary period so a future differently-periodic angle parameter needs no new formula, only a
    % new componentModelSpec.m .wrapParamPeriod entry.
    v = v - period*ceil((v - period/2)/period);
end
