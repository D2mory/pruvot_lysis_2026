function [counts, centers] = hexbin_plot(x, y, ncols, rngs)
% HEXBIN_PLOT  hexagonal binning and plot
%   hexbin_plot(x,y,ncols,rngs) builds a hex grid with approximately ncols
%   columns across the axis ranges rngs and plots the count per hexagon.
% Example:
%   x = randn(5000,1);
%   y = randn(5000,1);
%   hexbin_plot(x,y,40,[0 1 0 1]);

x = x(:);
y = y(:);

% bounding box
xmin=rngs(1);
xmax=rngs(2);
ymin=rngs(3);
ymax=rngs(4);

% choose hex side length s from desired number of columns
% horizontal spacing between centers = 1.5 * s
s = (xmax - xmin) / (1.5 * max(ncols-1,1));

% compute grid size
horiz = 1.5 * s;               % x spacing
vert  = sqrt(3) * s;           % y spacing

% determine number of columns and rows needed to cover bbox
ncol = ceil((xmax - xmin) / horiz) + 2;
nrow = ceil((ymax - ymin) / vert) + 2;

% generate hex centers (flat-top hexagons)
cols = (0:ncol-1);
rows = (0:nrow-1);
[C, R] = meshgrid(cols, rows);
xc = xmin + C(:) * horiz;
yc = ymin + R(:) * vert + (mod(C(:), 2) * (vert/2)); % odd columns offset

% trim centers outside bounding box
inside = (xc >= xmin - horiz) & (xc <= xmax + horiz) & ...
         (yc >= ymin - vert)  & (yc <= ymax + vert);
xc = xc(inside);
yc = yc(inside);
centers = [xc, yc];
ncent = size(centers,1);

% assign each point to the nearest center
% (Voronoi / nearest-center tessellation)
D = pdist2([x y], centers);
[~, idx] = min(D, [], 2);

% count points per center
counts = accumarray(idx, 1, [ncent, 1]);

% build one hexagon in local coordinates for plotting
theta = (0:5)*(2*pi/6);        % 6 vertices
% vertices coordinates = s * [cos(theta); sin(theta)]
hx = s*cos(theta);
hy = s*sin(theta);

% construct colormap
col=[.9 .2 .2];
col(2,:)=.7+.3*col;
itp=linspace(1,0,100).^3;
cmap=ones(100,1)*col(1,:)+itp'*(col(2,:)-col(1,:));
colormap(cmap)

% draw hexagons; only draw those with >0 points
for k = 1:ncent
    cx = centers(k,1);
    cy = centers(k,2);
    xv = cx + hx;
    yv = cy + hy;
    if counts(k) > 0
        patch(xv, yv, counts(k), 'EdgeColor', 0.8*[1 1 1]);
    else
        patch(xv, yv, [1 1 1], 'FaceAlpha', 0, 'EdgeColor', 0.8*[1 1 1]);
    end
end

xlim([xmin, xmax]);
ylim([ymin, ymax]);
set(gca,'box','on','layer','top')

end