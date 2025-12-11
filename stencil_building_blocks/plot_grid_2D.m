function plot_grid_2D(GRID,plot_centroids,varargin)
hold on;
if (GRID.twod)
    for n = 1,GRID.nblocks
        plot(GRID.gblock(n).x(:,:,1),GRID.gblock(n).y(:,:,1),varargin{:})
        plot(GRID.gblock(n).x(:,:,1).',GRID.gblock(n).y(:,:,1).',varargin{:})
    end
    if (plot_centroids)
        for n = 1,GRID.nblocks
            plot(GRID.gblock(n).xc,GRID.gblock(n).yc,'k.')
        end
    end
else
    error("That's not a 2D grid!");
end