function drawModelImage(ax, x, y, Z, cLim, cmap)
% DRAWMODELIMAGE  Draws a 2D scalar field (a fitted model evaluated on some grid -- isochromat- or
% voxel-resolution, either works) as an image into a GIVEN axes -- same "caller owns the axes"
% convention as plotVesselPatch.m, sharing that panel's own CLim/no-ticks styling (grayscale colormap by
% default) so a model panel reads as directly comparable to the patch panel it came from. Promoted out
% of fitVesselDiag.m's
% own panel-2 inline isochromat-image code (2026-08-18), generalized to also serve
% plotGaussianFitPanels.m's own downsampled/voxel-resolution model panel -- same drawing mechanics
% either way, only the input grid's resolution differs.
%
% INPUT
%   ax   : target axes (an existing axes handle).
%   x, y : 1D coordinate vectors (mm), matching imagesc's own axis-extent convention -- e.g.
%          XIso(1,:)/YIso(:,1) for an isochromat-resolution grid, or plotVesselPatch.m's own
%          toX(1:nx)/toY(1:ny) for a voxel-resolution one.
%   Z    : [numel(y) x numel(x)] scalar field to display.
%   cLim : [min max] -- pass plotVesselPatch.m's own returned cLim so this panel shares the SAME color
%          scale as the patch panel it's being compared against, rather than autoscaling independently.
%   cmap : [k x 3] (or []) colormap -- default gray(256) (this file's own original, only behavior).
%          Pass something else (e.g. divergingColormap.m) for a field that isn't itself an intensity on
%          the same scale as the patch it's drawn beside -- e.g. plotGaussianFitPanels.m's own residual
%          panel, a signed quantity that should never share the patch's plain grayscale ramp.
%
%   [toX, toY, cLim] = plotVesselPatch(ax1, im, voxSz2);
%   drawModelImage(ax2, XIso(1,:), YIso(:,1), Zfit, cLim);                       % grayscale (default)
%   drawModelImage(ax3, toX(1:nx), toY(1:ny), residual, cLimResidual, divergingColormap());  % diverging

    if nargin<6 || isempty(cmap); cmap = gray(256); end
    imagesc(ax, x, y, Z); axis(ax,'image'); colormap(ax,cmap);
    set(ax, 'CLim',cLim, 'XTick',[], 'YTick',[]);
end
