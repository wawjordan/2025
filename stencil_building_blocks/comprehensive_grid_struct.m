function GRID = comprehensive_grid_struct(file_name)
GRID_tmp = read_grd_file_to_struct(file_name);
GRID = grid_struct_check(GRID_tmp,true);
end