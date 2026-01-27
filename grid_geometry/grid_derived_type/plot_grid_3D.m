function plot_grid_3D(GRID,twod,plot_centroids,varargin)
hold on;
if twod
    for n = 1:GRID.nblocks
        plot_2D_grid_block(GRID.gblock(n),plot_centroids,varargin{:});
    end
else
    for n = 1:GRID.nblocks
        plot_3D_grid_block(GRID.gblock(n),plot_centroids,varargin{:});
    end
end
end

function plot_2D_grid_block(gblock,plot_centroids,varargin)
hold on;
plot( gblock.x(:,:,1)  , gblock.y(:,:,1),  varargin{:} )
plot( gblock.x(:,:,1).', gblock.y(:,:,1).',varargin{:} )
if (plot_centroids)
    p = plot(nan,nan,varargin{:}); color = get(p,'Color');
    plot( squeeze( gblock.grid_vars.cell_c(1,:,:,:) ), ...
          squeeze( gblock.grid_vars.cell_c(2,:,:,:) ),'.','Color',color)
end
end

function plot_3D_grid_block(gblock,plot_centroids,varargin)
hold on;
plot3Dgrid( gblock.x, gblock.y, gblock.z, varargin{:} );
if (plot_centroids)
    p = plot(nan,nan,varargin{:}); color = get(p,'Color');
    plot3( squeeze( gblock.grid_vars.cell_c(1,:,:,:) ), ...
           squeeze( gblock.grid_vars.cell_c(2,:,:,:) ), ...
           squeeze( gblock.grid_vars.cell_c(3,:,:,:) ),'.','Color',color)
end
end

function plot3Dgrid(X,Y,Z,varargin)
hold on;
sz = size(X);
for k = 1:sz(3)
    plot3( X(:,:,k),  Y(:,:,k),  Z(:,:,k),  varargin{:} )
end
for k = 1:sz(3)
    plot3( X(:,:,k).',Y(:,:,k).',Z(:,:,k).',varargin{:} )
end
for k = 1:sz(1)
    plot3( squeeze(X(k,:,:)).', ...
           squeeze(Y(k,:,:)).', ...
           squeeze(Z(k,:,:)).',varargin{:})
end

end