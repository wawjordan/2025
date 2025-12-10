classdef connected_bc_info_t
    properties
        my_block_id(1,1)   int32 = -1
        my_face_id(1,1)    int32 = -1
        my_bnd_min(3,1)    int32 = -1
        my_bnd_max(3,1)    int32 = -1
        nbor_block_id(1,1) int32 = -1
        nbor_face_id(1,1)  int32 = -1
        nbor_bnd_min(3,1)  int32 = -1
        nbor_bnd_max(3,1)  int32 = -1
    end
end