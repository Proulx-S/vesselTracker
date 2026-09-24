function [hMeas, hPred] = drawRadialProfileMeasuredPredicted(ax, teff, measured, predicted, color)
% DRAWRADIALPROFILEMEASUREDPREDICTED  Draws one fit's per-voxel measured/predicted intensity, each
% against its own effective distance from center (see effectiveDistanceFromFit.m), as two scatter
% series into a GIVEN axes -- same "caller owns the axes" convention as plotVesselPatch.m. Promoted out
% of fitVesselDiag.m's own panel-3 inline code (2026-08-18), split from
% drawRadialProfileModelContinuous.m so a caller can show either series independently (e.g.
% plotGaussianFitPanels.m's own separate measured/predicted vs. model-construction panels) or both
% together on the SAME axes (fitVesselDiag.m's own panel 3, which draws BOTH this and the continuous
% model curve on one axes).
%
% INPUT
%   ax        : target axes.
%   teff      : [n x 1] effective distance from center (mm), one per fit voxel.
%   measured  : [n x 1] that voxel's own raw intensity.
%   predicted : [n x 1] (or []) that voxel's own isochromat-averaged predicted intensity (SAME order as
%               teff) -- pass [] to skip the predicted series entirely (e.g. plotGaussianFitPanels.m's
%               own xd panel, which shows measured data only).
%   color     : RGB triplet, shared by both series (measured filled/outlined in black, predicted open).
%
% OUTPUT
%   hMeas, hPred : scatter handles, for a caller building its own legend. hPred is gobjects(0) if
%                  predicted was empty/omitted.
%
%   [hMeas,hPred] = drawRadialProfileMeasuredPredicted(ax, teff, measured, predicted, [1 0.85 0.1]);
%   hMeas         = drawRadialProfileMeasuredPredicted(ax, teff, measured, [], [1 1 1]);   % measured only

    hold(ax,'on'); set(ax,'Color','k','XColor','w','YColor','w');
    hMeas = scatter(ax, teff, measured, 40, color, 'filled', 'MarkerEdgeColor','k');
    hPred = gobjects(0);
    if nargin>=4 && ~isempty(predicted)
        hPred = scatter(ax, teff, predicted, 18, color);
    end
end
