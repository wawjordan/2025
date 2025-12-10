classdef grid_block_lite

    properties

        % imax, jmax, and kmax are the number of nodes in the block
        imax    (1,1) {mustBeInteger}
        jmax    (1,1) {mustBeInteger}
        kmax    (1,1) {mustBeInteger}

        % i_cells, j_cells, and k_cells are the number of cells in the block
        i_cells (1,1) {mustBeInteger}
        j_cells (1,1) {mustBeInteger}
        k_cells (1,1) {mustBeInteger}

        Nnodes  (1,3) {mustBeInteger}
        Ncells  (1,3) {mustBeInteger}

        % sized nblocks, set
        total_cells  (1,1) {mustBeInteger}
        total_volume (1,1) {mustBeNumeric}

        % parent block id
        block_id     (1,1) {mustBeInteger}
        % x,y,z nodal coordinates for each block, i.e., x(i,j,k)
        x (:,:,:) {mustBeNumeric}
        y (:,:,:) {mustBeNumeric}
        z (:,:,:) {mustBeNumeric}

        allocated (1,1) {mustBeNumericOrLogical} = false

    end

    methods
        function this = allocate_grid_block( this, imax, jmax, kmax )
            this.imax = imax;
            this.jmax = jmax;
            this.kmax = kmax;
            this.i_cells = imax-1;
            this.j_cells = jmax-1;
            this.k_cells = kmax-1;
            this.Nnodes = [imax;jmax;kmax];
            this.Ncells = [imax-1;jmax-1;kmax-1];
            this.x = zeros( imax, jmax, kmax );
            this.y = zeros( imax, jmax, kmax );
            this.z = zeros( imax, jmax, kmax );
            this.allocated = true;
        end

        function sub_grid = subset_grid_block(this,lo,hi)
            sub_grid = grid_block_lite();
            sub_grid.imax = hi(1) - lo(1) + 1;
            sub_grid.jmax = hi(2) - lo(2) + 1;
            sub_grid.kmax = hi(3) - lo(3) + 1;
            sub_grid.i_cells = sub_grid.imax-1;
            sub_grid.j_cells = sub_grid.jmax-1;
            sub_grid.k_cells = sub_grid.kmax-1;
            sub_grid.Nnodes = [sub_grid.imax;sub_grid.jmax;sub_grid.kmax];
            sub_grid.Ncells = [sub_grid.imax-1;sub_grid.jmax-1;sub_grid.kmax-1];
            sub_grid.x  = this.x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
            sub_grid.y  = this.y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
            sub_grid.z  = this.z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
            sub_grid.allocated = true;
        end

    end % methods

    methods (Static)
        function [Xtmp,Ytmp,Ztmp] = copy_gblock_nodes(gblock,cell_idx,nSkip)
            Xtmp = zeros(nSkip(1)+1,nSkip(2)+1,nSkip(3)+1);
            Ytmp = zeros(nSkip(1)+1,nSkip(2)+1,nSkip(3)+1);
            Ztmp = zeros(nSkip(1)+1,nSkip(2)+1,nSkip(3)+1);

            nCellNodes = nSkip + 1;
            nPoints    = prod(nCellNodes);
            for j = 1:nPoints
                offset = grid_block.global2local( j, nCellNodes );
                itmp   = cell_idx + offset - 1;
                Xtmp(offset(1),offset(2),offset(3)) = gblock.x(itmp(1),itmp(2),itmp(3));
                Ytmp(offset(1),offset(2),offset(3)) = gblock.y(itmp(1),itmp(2),itmp(3));
                Ztmp(offset(1),offset(2),offset(3)) = gblock.z(itmp(1),itmp(2),itmp(3));
            end
        end
        function iSub = global2local(iG,nSub)
            nDims = numel(nSub);
            iSub = zeros(1,nDims);
            if (nDims==1)
                iSub(1) = iG;
                return
            end
            p = prod(nSub);
            iGtmp = iG;
            for i = nDims:-1:1
                p = fix( p/nSub(i) );
                iTmp = mod(iGtmp-1,p)+1;
                iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
                iGtmp = iTmp;
            end
        end
    end
end