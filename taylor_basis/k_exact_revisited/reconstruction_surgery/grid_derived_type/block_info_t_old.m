classdef block_info_t_old
    properties
        block_id(1,1)   int32 = 0
        Ncells(3,1)     int32 = 0
        connected_bc_list(:,1) connected_bc_info_t_old
    end
end