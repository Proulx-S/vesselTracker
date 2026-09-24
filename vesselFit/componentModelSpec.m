function spec = componentModelSpec(name)
% COMPONENTMODELSPEC  The MODEL REGISTRY for fitComponentFit.m -- everything the solver needs to know
% about one kind of additive model COMPONENT, returned as plain data.
%
% A fit is a LIST OF COMPONENTS whose predictions are SUMMED:
%     pred = sum_c  modelFun_c(X, Y, toModelParams_c(vals_c){:})
% Every component is one entry of a fitComponentFit.m components array; this function turns a model
% NAME into the spec that entry needs. Adding a new model type (e.g. 'parab') is a new case here plus
% a formula file -- fitComponentFit.m itself never branches on model type at all.
%
% WHY 'gaussian' HAS NO 'b' -- unlike fitVessel.m's own 7-parameter 'gaussian' method, the component
% gaussian is SIX parameters (a,x0,y0,radius,aspectRatio,theta). The additive background is its own
% separate 'background' component, fit JOINTLY with the peaks rather than carried inside each of them.
% Two reasons this is the right split, not just a rearrangement:
%   1. With K peaks sharing one patch there is physically ONE background, not K of them. Carrying a
%      per-peak `b` forces either K redundant baselines (over-parameterized, and they trade off against
%      each other) or an out-of-band "which peak's b is the real one" rule. fitMultiVessel.m already
%      had to solve this with a bespoke one-shared-b-plus-b-union-bound arrangement; making background a
%      component makes that structural instead of special-cased.
%   2. The background then gets the SAME per-parameter granularity surface every peak parameter has, for
%      free -- 'perRun' or 'perFrame' baseline drift is expressible without any new machinery, where
%      fitVessel.m's 7-parameter gaussian could only ever say "b, like every other parameter".
% gaussianModel.m ITSELF is untouched (still takes a literal `b` positionally) -- 'gaussian's own
% toModelParams simply passes b=0, and the 'background' component supplies the real one. See
% fitComponentFit.m's own LINEARITY note for why summing them is exact.
%
% PARAMETERIZATION -- 'gaussian' optimizes (radius,aspectRatio), never (sx,sy). There is no
% parameterization CHOICE anywhere in this project (see fitVessel.m's own PARAMETERIZATION note);
% radius=sqrt(sx*sy) and aspectRatio=sx/sy are converted to gaussianModel.m's own literal sx,sy inside
% toModelParams, once per evaluation, as plain unnamed locals -- never stored, never reported.
%
% NO NO-ARG CALL, deliberately -- a pure name->data lookup with no opts struct, same documented
% exemption from the self-populating-default-opts convention as gaussianModel.m/buildIsochromatGrid.m/
% additiveSeedBound.m. The opts-bearing entry point is fitComponentFit.m, which does have one.
%
% INPUT
%   name : 'gaussian' | 'background'
%
% OUTPUT
%   spec : struct --
%     .name             the name, echoed back (provenance).
%     .paramNames       cellstr, THE canonical order for this model. Every per-parameter struct a
%                       caller supplies (.mode/.fixed/.start/.lower/.upper) is keyed by these names,
%                       and fitComponentFit.m packs its free vector in this order.
%     .modelFun         @(X,Y,args...) -> value on the isochromat grid.
%     .toModelParams    @(vals) -> cell of modelFun's own positional args, from a cell in .paramNames
%                       order. This is the ONLY place a parameterization/argument-order difference
%                       between the optimizer's view and the formula's view is allowed to live.
%     .perPeak          logical -- does ONE entry of this model expand into one component PER SEED PEAK
%                       (true, 'gaussian'), or is it structurally SINGLETON regardless of peak count
%                       (false, 'background')? fitVesselPatchTimeSeries.m reads this to decide two
%                       things: whether an entry participates in the peak-CLAIMING accounting that
%                       makes "every remaining peak" work, and whether its seed selects peaks to model
%                       or merely names where to read a starting value from. A singleton model that
%                       claimed peaks would starve the per-peak entries (or, if it ran last, find every
%                       peak already claimed and resolve to nothing) -- so this is load-bearing, not
%                       descriptive.
%     .shapeParamNames  the two shape parameters, or {} for a model with no ellipse geometry. When
%                       non-empty the model is expected to also have a 'theta'.
%     .canonicalizeShapeFun / .shapeSwapIdxFun / .canonicalizeShapeBoundsFun / .shapeBoundsReciprocal
%                       post-fit degeneracy reduction -- IDENTICAL contract and semantics to
%                       fitVessel.m's own gspec fields of the same names (see gaussianMethodSpec there
%                       for the full derivation of each). Empty/absent for a shapeless model.
%     .wrapParamPeriod  struct, paramName -> scalar period (radians), for any parameter that is
%                       genuinely PERIODIC rather than a plain unbounded real -- absent/no field for
%                       every other parameter. fitComponentFit.m consults this ONLY when expanding a
%                       mode='perRunPoly' parameter's fitted within-run polynomial back to a per-frame
%                       trajectory (see that file's own mode='perRunPoly' doc): the raw polynomial value
%                       can legitimately wander outside one period as the trajectory drifts, so the
%                       EXPANDED per-frame report is wrapped into (-period/2, period/2]. 'gaussian'.theta
%                       is period pi (NOT 2*pi): canonicalizeRadiusAspectShape's own aspectRatio<->1/
%                       aspectRatio swap already treats theta and theta+pi as the identical ellipse, and
%                       this reuses that EXACT SAME wrap formula (see canonicalizeRadiusAspectShape
%                       below), just generalized over an arbitrary period so a future genuinely-2*pi
%                       angle parameter needs no new formula, only a new struct entry here.
%
%   spec = componentModelSpec('gaussian');
%
% See also fitComponentFit, gaussianModel, canonicalizeEllipseShape, fitVessel.

    switch name
        case 'gaussian'
            spec.name           = 'gaussian';
            spec.paramNames     = {'a','x0','y0','radius','aspectRatio','theta'};
            spec.perPeak        = true;    % one component per seed peak
            spec.modelFun       = @gaussianModel;
            spec.toModelParams  = @gaussianComponentToModelParams;
            spec.shapeParamNames = {'radius','aspectRatio'};
            spec.canonicalizeShapeFun = @canonicalizeRadiusAspectShape;
            % aspectRatio<1 is the reciprocal-swap test -- the SAME predicate
            % canonicalizeRadiusAspectShape makes internally. s1=radius is swap-invariant and unused.
            spec.shapeSwapIdxFun = @(s1,s2) s2 < 1;
            % IDENTITY, not a swap: radius and aspectRatio are not interchangeable quantities (a length
            % against a dimensionless ratio), so swapping their BOUND fields would put radius-scale
            % numbers under aspectRatio's row -- a real, confirmed bug in fitVessel.m's own history
            % (2026-08-17). Once aspectRatio's bound is reciprocal-symmetric around 1, that interval is
            % invariant under the reciprocal anyway, so there is nothing to update.
            spec.canonicalizeShapeBoundsFun = @(v1,v2,swapIdx) deal(v1,v2);
            spec.shapeBoundsReciprocal = true;
            spec.wrapParamPeriod = struct('theta', pi);
        case 'background'
            % A single additive constant over the whole patch. Its "model" is deliberately a real
            % function of the grid rather than a scalar special case, so fitComponentFit.m can treat
            % every component through one uniform evaluate-and-sum path with no isbackground branch.
            spec.name           = 'background';
            spec.paramNames     = {'b'};
            spec.perPeak        = false;   % SINGLETON -- one background whatever the peak count
            spec.modelFun       = @(X, ~, b) b * ones(size(X));
            spec.toModelParams  = @(vals) vals;   % identity -- one parameter, already in arg order
            spec.shapeParamNames = {};
            spec.canonicalizeShapeFun = [];
            spec.shapeSwapIdxFun = [];
            spec.canonicalizeShapeBoundsFun = [];
            spec.shapeBoundsReciprocal = false;
            spec.wrapParamPeriod = struct();
        otherwise
            error('componentModelSpec:unknownModel', ...
                ['unknown component model ''%s'' -- registered models are ''gaussian'' and ' ...
                 '''background''. Add a case here (plus its formula file) to register a new one; ' ...
                 'fitComponentFit.m itself needs no change.'], char(string(name)));
    end
end

% ---------------------------------------------------------------------------
function modelVals = gaussianComponentToModelParams(vals)
    % vals is a 1x6 cell in paramNames order (a,x0,y0,radius,aspectRatio,theta) -> gaussianModel.m's own
    % literal (a,x0,y0,sx,sy,theta,b) argument order, with b FORCED TO ZERO: the additive background is
    % a separate component in this framework (see file header), so a gaussian peak here contributes its
    % peak shape ONLY. Passing the real background in here instead would double-count it once per peak.
    % sx=radius*sqrt(aspectRatio), sy=radius/sqrt(aspectRatio) -- the exact inverse of
    % radius=sqrt(sx*sy), aspectRatio=sx/sy.
    sqrtAR = sqrt(vals{5});
    sx = vals{4} * sqrtAR;
    sy = vals{4} / sqrtAR;
    modelVals = {vals{1}, vals{2}, vals{3}, sx, sy, vals{6}, 0};
end

% ---------------------------------------------------------------------------
function [radius, aspectRatio, theta] = canonicalizeRadiusAspectShape(radius, aspectRatio, theta)
    % Duplicated from fitVessel.m's own canonicalizeRadiusAspectShape (this codebase's established
    % small-helper-duplication convention -- fitMultiVessel.m carries the same copy). radius is
    % swap-invariant (sqrt(sx*sy)=sqrt(sy*sx)), so only aspectRatio needs the analog of
    % canonicalizeEllipseShape.m's s1<->s2 exchange -- and exchanging a ratio with a length is
    % meaningless, so the correct analog is a RECIPROCAL: aspectRatio<1 means sy was the larger, and the
    % identical ellipse is equally described by 1/aspectRatio with theta a further quarter-turn round.
    swapIdx = aspectRatio < 1;
    aspectRatio(swapIdx) = 1 ./ aspectRatio(swapIdx);
    theta(swapIdx) = theta(swapIdx) + pi/2;
    theta = theta - pi*ceil((theta - pi/2)/pi);   % same period-pi wrap as canonicalizeEllipseShape.m
end
