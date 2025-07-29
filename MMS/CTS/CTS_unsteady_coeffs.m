%% CTS Unsteady NS: 2D coefficients
clc; clear;
a = zeros(22,5);

lref = 1.0;
dest='C:\Users\Will\Documents\MATLAB\VT_Research\2022_SUMMER\anisotropic\test_geometry\curvilinear_mms\';
filename = 'mms.config';

%% rho
a( 1,1) =  1.0;
a( 2,1) =  0.15;
a( 3,1) =  2.0;
a( 4,1) =  1/3;
a( 5,1) = -0.1;
a( 6,1) =  1.0;
a( 7,1) = -0.2;
a( 8,1) =  0.0;
a( 9,1) =  0.0;
a(10,1) =  0.0;
a(11,1) = -0.05;
a(12,1) =  1.0;
a(13,1) =  0.25;
a(14,1) =  0.0;
a(15,1) =  0.0;
a(16,1) =  0.0;
a(17,1) =  0.0;
a(18,1) =  0.0;
a(19,1) =  0.0;
a(20,1) =  0.1;
a(21,1) =  2.0;
a(22,1) =  1.0;

%% u
a( 1,2) =  40.0;
a( 2,2) =  2.5;
a( 3,2) =  2/3;
a( 4,2) =  0.5;
a( 5,2) = -1.5;
a( 6,2) =  1.5;
a( 7,2) = -1/8;
a( 8,2) =  0.0;
a( 9,2) =  0.0;
a(10,2) =  0.0;
a(11,2) = -3.0;
a(12,2) =  1.5;
a(13,2) =  1/8;
a(14,2) =  0.0;
a(15,2) =  0.0;
a(16,2) =  0.0;
a(17,2) =  0.0;
a(18,2) =  0.0;
a(19,2) =  0.0;
a(20,2) =  2.0;
a(21,2) =  2/3;
a(22,2) =  2.0;

%% v
a( 1,3) =  40.0;
a( 2,3) = -3.8;
a( 3,3) =  5/3;
a( 4,3) = -0.5;
a( 5,3) =  2.0;
a( 6,3) =  3.0;
a( 7,3) =  0.0;
a( 8,3) =  0.0;
a( 9,3) =  0.0;
a(10,3) =  0.0;
a(11,3) =  3.0;
a(12,3) =  1.5;
a(13,3) =  1/8;
a(14,3) =  0.0;
a(15,3) =  0.0;
a(16,3) =  0.0;
a(17,3) =  0.0;
a(18,3) =  0.0;
a(19,3) =  0.0;
a(20,3) =  1.5;
a(21,3) =  1.25;
a(22,3) =  2.0;

%% w
a( :,4) = 0.0;

%% p
a( 1,5) =  1.0E5;
a( 2,5) = -0.1E5;
a( 3,5) =  0.5;
a( 4,5) =  0.5;
a( 5,5) =  0.15E5;
a( 6,5) =  1.0;
a( 7,5) = -0.5;
a( 8,5) =  0.0;
a( 9,5) =  0.0;
a(10,5) =  0.0;
a(11,5) = -0.1E5;
a(12,5) =  0.75;
a(13,5) = -0.25;
a(14,5) =  0.0;
a(15,5) =  0.0;
a(16,5) =  0.0;
a(17,5) =  0.0;
a(18,5) =  0.0;
a(19,5) =  0.0;
a(20,5) =  0.1E5;
a(21,5) =  1.0;
a(22,5) =  0.25;


write_config(a,lref,fullfile(dest,filename));

function write_config(A,lref,filename)

fltfmt = '%+#23.16E  ';

max_mms_num  = 1;
file_nlines  = 5;
file_ninputs = 23;

mms_num = 1;

fid = fopen(filename,'w');

fprintf(fid,'  %s\n','laminar_ns');
fprintf(fid,'  %d  %d  %d\n',max_mms_num, file_nlines, file_ninputs);
fprintf(fid,'  %d\n',mms_num);

for i = 1:file_nlines
    fprintf(fid,[repmat(fltfmt,1,file_ninputs-1),'%f\n'],[A(:,i);lref]);
end

fclose(fid);
end