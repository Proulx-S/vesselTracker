function [hLine, hIso] = drawRadialProfileModelContinuous(ax, tCurve, curveVals, teffIso, valsIso, color)
% DRAWRADIALPROFILEMODELCONTINUOUS  Draws a fit's own radial profile as a CONTINUOUS function
% (diag.modelFun's 1D mode, see gaussianModel.m/parabolicModel.m's own 1D MODE docs), optionally
% overlaid with small markers at individual ISOCHROMAT sample points (the sub-voxel grid the fit was
% actually evaluated/averaged on before reduction to voxel resolution -- see buildIsochromatGrid.m) --
% into a GIVEN axes, same "caller owns the axes" convention as plotVesselPatch.m. Promoted out of
% fitVesselDiag.m's own panel-3 inline code (2026-08-18), split from
% drawRadialProfileMeasuredPredicted.m -- see that file's own header for why.
%
% INPUT
%   ax        : target axes.
%   tCurve    : [1 x m] (or []) effective-distance sample points for the continuous curve (e.g.
%               linspace(0,teffMax,200)) -- pass [] to skip the curve entirely (e.g.
%               plotGaussianFitPanels.m's own xe panel, which shows isochromat markers only, no line).
%   curveVals : [1 x m] diag.modelFun(tCurve,[],modelVals{:}) -- the model evaluated at tCurve. Ignored
%               if tCurve is empty.
%   teffIso   : [k x 1] (or []) effective distance of individual isochromat samples -- pass [] to
%               skip the marker overlay entirely (e.g. fitVesselDiag.m's own panel 3, which shows only
%               the continuous curve, no per-isochromat markers).
%   valsIso   : [k x 1] the model's own value at each teffIso sample (SAME order). Ignored if teffIso is
%               empty.
%   color     : RGB triplet, shared by the curve and the isochromat markers.
%
% OUTPUT
%   hLine : the continuous-curve line handle, or gobjects(0) if tCurve was empty/omitted.
%   hIso  : the isochromat-marker scatter handle, or gobjects(0) if teffIso was empty/omitted.
%
%   [hLine,~]     = drawRadialProfileModelContinuous(ax, tCurve, curveVals, [], [], [1 0.85 0.1]);          % curve only
%   [hLine,hIso]  = drawRadialProfileModelContinuous(ax, tCurve, curveVals, teffIso, valsIso, [1 0.85 0.1]); % + isochromat markers
%   [~,hIso]      = drawRadialProfileModelContinuous(ax, [], [], teffIso, valsIso, [1 1 1]);                % isochromat markers only

    hold(ax,'on'); set(ax,'Color','k','XColor','w','YColor','w');
    hLine = gobjects(0);
    if ~isempty(tCurve)
        hLine = plot(ax, tCurve, curveVals, '-', 'Color',color, 'LineWidth',1.5);
    end
    hIso = gobjects(0);
    if nargin>=5 && ~isempty(teffIso)
        hIso = scatter(ax, teffIso, valsIso, 6, color, 'filled', 'MarkerFaceAlpha',0.35);
    end
end
