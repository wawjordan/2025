function DATA = tec2mat_structured_condensed(file_name)

ZONES = tec2mat_MOD2_local(file_name);
DATA = struct();
for i = 1:numel(ZONES)
    DATA_TMP = get_structured_from_sensei_tec2mat_local(ZONES,i);
    DATA(i).ZONE = make_named_struct_local(DATA_TMP);
end

end

function DATA2 = make_named_struct_local(DATA)
N = numel(DATA);

inputs = cell(1,2*N);

cnt = 0;
for i = 1:N
    cnt = cnt + 1;
    % inputs{cnt} = replace(DATA(i).name,{' ','(',')','-'},'_');
    inputs{cnt} = matlab.lang.makeValidName(DATA(i).name);
    cnt = cnt + 1;
    inputs{cnt} = DATA(i).data;
end
    
DATA2 = struct(inputs{:});
end

function DATA = get_structured_from_sensei_tec2mat_local(ZONES,zone_num)

ZONE = ZONES(zone_num);
nvars = length(ZONE.VARLOC);
DATA = struct();
cell_var_cnt = 0;
node_var_cnt = 0;

% imax = ZONE.I;
if isfield(ZONE,'J')
    jmax = ZONE.J;
else
    jmax = 1;
end

if isfield(ZONE,'K')
    kmax = ZONE.K;
else
    kmax = 1;
end

o = 0;
for n = 1:nvars
    DATA(n+o).name = ZONE.ZONEVARlist{n};
    if ZONE.VARLOC(n) == 0 % nodal
        node_var_cnt = node_var_cnt + 1;
        if jmax > 1
            if kmax > 1
                DATA(n+o).data = ZONE.data_nodal(:,:,:,node_var_cnt);
            else
                DATA(n+o).data = ZONE.data_nodal(:,:,node_var_cnt);
            end
        else
            DATA(n+o).data = ZONE.data_nodal(:,node_var_cnt);
        end
        % add cell-centered coordinates
        if any(strcmp(DATA(n+o).name,{'X','Y','Z'}))
            hi = size(DATA(n+o).data,[1,2,3]);
            o = o + 1;
            DATA(n+o).name = [DATA(n+o-1).name,'C'];
            DATA(n+o).data = zeros(hi(1)-1,hi(2)-1,hi(3)-1);
            for k = 1:hi(3)-1
                for j = 1:hi(2)-1
                    for i = 1:hi(1)-1
                        DATA(n+o).data(i,j,k) = mean(DATA(n+o-1).data(i:i+1,j:j+1,k:k+1),"all");
                    end
                end
            end
        end         
    elseif ZONE.VARLOC(n) == 1 % cell-centered
        cell_var_cnt = cell_var_cnt + 1;
        if jmax > 2
            if kmax > 2
                DATA(n+o).data = ZONE.data_cell(:,:,:,cell_var_cnt);
            else
                DATA(n+o).data = ZONE.data_cell(:,:,cell_var_cnt);
            end
        else
            DATA(n+o).data = ZONE.data_cell(:,cell_var_cnt);
        end
    else
        error('VARLOC should either be 0 or 1')
    end
end

end

function [zone,VARlist,auxdata] = tec2mat_MOD2_local(fname,varargin)
%% load data file in Tecplot ASCII format
%% Author: Mach6
% Github repository: https://github.com/luohancfd/FluidDynamicTools/tree/master/Tecplot_Tools
% Method of calling:
%   [zone,VARlist] = tec2mat(fname,debug);
%   fname(in): name of data file
%   optional flags:
%          "debug" print zone information out (default: false)
%          "safe" check the existence of scientific notation like 100-100 (default: false)
%        For IJK ordered date only:
%          "nopassive" include passive variables (0) in the data field (default: true)
%          "passive"   not include passive variables in the data field
% -------------------------------------------------------------------------
% ChangeLog:
% v1: support IJK data without mixing data format (node with cell-centered)
%     support Finite Element grid
% v2: add a "safe" optional flag to handle number like 100-100 (100E-100)
%--------------------------------------------------------------------------
debug = false;
nopassive = true;
safe = false;
if nargin > 1
    for i = 1:length(varargin)
        if strcmp(varargin{i},'debug')
            debug = true;
        elseif strcmp(varargin{i},'nopassive')
            nopassive = true;
        elseif strcmp(varargin{i},'passive')
            nopassive = false;
        elseif strcmp(varargin{i},'safe')
            safe = true;
        else
            error('Unknown option: %s',varargin{i})
        end
    end
end
%% First test if the file has "strange" scientific notations like  2.3-100 (2.3E-100)
if ~safe
    IMY_FSCANF = false;
else
    filetext = fileread(fname);
    strange = regexp(filetext,'(?<!(?:(,|\[)))-?(\d+(\.\d+)?|\.\d+)([-+]\d+)','match','once');
    if isempty(strange)
        IMY_FSCANF = false;
    else
        IMY_FSCANF = true;
    end
end
%% load tecplot ASCII file and write result to zone and VARList
fid = fopen(fname,'r');
NVAR = 0;
%% read file untill we get a line of variables name
IVAR = false;
varline = '';
while ~feof(fid)
    pos = ftell(fid);  %line above ZONE
    line = upper( strtrim( fgetl(fid) ) );
    if ~isempty(line)
        if (~strcmp(line(1),'#'))
            % if contains(line,'VARIABLES')
            if ~isempty(regexp(line,'VARIABLES\s*=','match'))
                varnames = regexp(line,'VARIABLES\s*=\s*(.*)','tokens');
                varline = [varline, varnames{1}{1}];
                IVAR = true;
            elseif (contains(line,'ZONE') || contains(line,'DATASETAUXDATA'))
                fseek(fid,pos,-1); %jump out at position above zone
                break
            elseif (IVAR)
                varline = [varline, line];
            end
        end
    end
end

if debug
    fprintf('===========================================\n');
    fprintf('Read file: %s\n',fname);
end
if IVAR
    varline0 = varline;
    % replace all delimiter by a special symbol
    varline = strrep(varline,'''','$#');
    varline = strrep(varline,'"','$#');
    varline = strrep(varline,',','$#');
    VARlist0 =strsplit(varline,'$#');
    VARlist = {};
    for i = 1:length(VARlist0)
        if ~isempty(strtrim(VARlist0{i}))
            VARlist = [VARlist, VARlist0{i}];
        end
    end
    
    NVAR = length(VARlist);
    if debug
        fprintf('The line defining variables:\n%s\n',varline0);
        fprintf('The variables we get:\n');
        fprintf('%s, ',VARlist{:});
        fprintf('\n');
    end
end

nzone = 0;
zone = [];
auxdata = {};
if ~IVAR
    error("there is no variable list")
else
    while ~feof(fid)
        line = upper( strtrim( fgetl(fid) ) );
        if ~isempty(line)
            if strcmp(line(1),'#')
                continue
            end
        end
        if contains(line,'DATASETAUXDATA')
            if strcmpi(line(1:14),'DATASETAUXDATA')
                line=strtrim(line(15:end));
                temp = strsplit(line,'=');
                temp(2)=strtrim(temp(2));
                temp(2)=strrep(temp(2),'"','');
                temp(2)=strrep(temp(2),'''','');
                auxdata = [auxdata;{temp(1),temp(2)}];
            end
        end
            
        if contains(line,'ZONE')
            nzone = nzone + 1;
                           
            zoneline = regexp(line,'ZONE(.*)','tokens');
            zoneline = strtrim(zoneline{1}{1});
            
            while ~feof(fid)
                pos = ftell(fid);  %position above data
                line = upper( strtrim( fgetl( fid ) ) );
                if ~isempty(line)
                    if strcmp(line(1),'#')
                        continue
                    end
                end
                if contains(line,'=')
                    zoneline = [zoneline,' ',line];
                else
                    % arrive at data
                    fseek(fid,pos,-1); %jump out at position above zone
                    break
                end
            end
            
            zoneline0 = zoneline;
            % preprocess zone property, replace equal sign in quotes with a
            % special symbol
            quoted =  regexp(zoneline,'"(.+?)"','tokenExtents');
            quoteds = regexp(zoneline,'''(.+?)''','tokenExtents');
            quoted = [quoted quoteds];
            
            zoneline = '';
            iendp = 1;
            for i = 1:length(quoted)
                istart = quoted{i}(1);
                iend = quoted{i}(2);
                if istart-1 >= iendp
                    zoneline = [zoneline,zoneline0(iendp:istart-1)];
                end
                zoneline = [zoneline,strrep(zoneline0(istart:iend),'=','$#')];
                if iend == length(zoneline0)
                    break
                end
                iendp = iend+1;
            end
            if iend ~= length(zoneline0)
                zoneline = [zoneline, zoneline0(iendp:end)];
            end
                
            % get zone properties
            [ZONEprop.name,ZONEprop.start,ZONEprop.end] = regexp(zoneline,'\w+(?=\s*=)\s*=','match','start','end');
            Nzoneprop = length(ZONEprop.name);
            ZONEprop.prop = cell(1,Nzoneprop);
            for i = 1:Nzoneprop-1
                ZONEprop.name{i} = strtrim(ZONEprop.name{i}(1:end-1)); %remove "=" sign
                ZONEprop.name{i} = upper(ZONEprop.name{i});  %upper case
                if strcmpi(ZONEprop.name{i},'Nodes')
                    ZONEprop.name{i} = 'N';
                end
                if strcmpi(ZONEprop.name{i},'Elements')
                    ZONEprop.name{i} = 'E';
                end
                
                ZONEprop.prop{i} = strtrim(zoneline(ZONEprop.end(i)+1:ZONEprop.start(i+1)-1));
                if strcmp(ZONEprop.prop{i}(end),',')
                    ZONEprop.prop{i} = ZONEprop.prop{i}(1:end-1);
                end
            end
            ZONEprop.name{Nzoneprop} = strtrim(ZONEprop.name{Nzoneprop}(1:end-1)); 
            ZONEprop.prop{Nzoneprop} = zoneline(ZONEprop.end(Nzoneprop)+1:end);
            
            
            % reformat properties
            Prop = [];
            for i = 1:Nzoneprop
                if isstrprop(ZONEprop.prop{i},'digit')
                    Prop.(ZONEprop.name{i}) = str2num(ZONEprop.prop{i});
                else
                    Prop.(ZONEprop.name{i}) = strtrim(strrep(ZONEprop.prop{i},'$#','='));
                end
            end
            clear ZONEprop;
            
            if ~isfield(Prop,'T')
                Prop.T = sprintf('ZONE %d',nzone);
            else
                if strcmp(Prop.T(1),'"') && strcmp(Prop.T(end),'"')
                    Prop.T = Prop.T(2:end-1);
                elseif strcmp(Prop.T(1),'''') && strcmp(Prop.T(end),'''')
                    Prop.T = Prop.T(2:end-1);
                end
            end
            Prop.title = Prop.T;
                
            localVARindex = 1:NVAR;
            if isfield(Prop,'PASSIVEVARLIST')
                list = regexp(Prop.PASSIVEVARLIST,'[(.*)\]','tokens');
                if ~isempty(list)
                    PASSIVEVARLIST = strl2list(list{1}{1});
                    Prop.PASSIVEVARLIST = PASSIVEVARLIST;
                    localVARindex = setdiff(localVARindex,Prop.PASSIVEVARLIST);
                    ZONEVARlist = VARlist(localVARindex);
                else
                    ZONEVARlist = VARlist;
                end
            else
                Prop.PASSIVEVARLIST =[];
                ZONEVARlist = VARlist;
            end
            localNVAR = length(ZONEVARlist);
            
            VARLOC = zeros(1,NVAR);  %0 for node, 1 for cell
            Prop.IVARLOC = false;
            if isfield(Prop,'VARLOCATION')
                varloc = regexp(Prop.VARLOCATION,'\((.*)\)','tokens');
                if ~isempty(varloc)
                    Prop.IVARLOC = true;
%                     varloc = strsplit(varloc{1}{1},{','});
                    varloc = split(varloc{1}{1},regexpPattern('(,)(?=(\[\S*\]))'));
                    nodalvar = []; cellvar = [];
                    for i = 1:length(varloc)
                        if contains(varloc{i},'NODAL')
                            nodalvar =  regexp(varloc{i},'\[(.*)\]','tokens');
                            nodalvar = strl2list(nodalvar{1}{1});
                        elseif contains(varloc{i},'CELLCENTERED')
                            cellvar =  regexp(varloc{i},'\[(.*)\]','tokens');
                            cellvar = strl2list(cellvar{1}{1});
                        end
                    end
                    
                    for i = 1:NVAR
                        if ismember(i,cellvar)
                            VARLOC(i) = 1;
                        end
                    end
                end
            end
            Prop.VARLOC =VARLOC; 

            if isfield(Prop,'I')
                % IJK ordered ZONE
                if isfield(Prop,'ZONETYPE')
                    if ~strcmpi(Prop.ZONETYPE,'ORDERED')
                        error(['ZONE ',num2str(nzone),' type is not ordered']);
                    end
                else
                    Prop.ZONETYPE = 'ORDERED';
                end
                idata = Prop.I;
                jdata = 1; kdata = 1;
                if isfield(Prop,'J')
                    jdata = Prop.J;
                end                
                if isfield(Prop,'K')
                    kdata = Prop.K;
                end
                if debug
                    fprintf('ZONE %d is a IJK ordered zone\n',nzone);
                    fprintf('ZONE %d Size: %d x %d x %d NVar: %d\n',nzone,idata,jdata,kdata,localNVAR);
                end
                ndata = idata*jdata*kdata;
                if ~isfield(Prop,'DATAPACKING')
                    Prop.DATAPACKING = 'POINT';
                end

                % assume all nodal at first
                data_nodal = [];
                data_cell  = [];

                ndata_total = ndata*localNVAR;

                Ncells(1) = idata-1;
                Ncells(2) = jdata;
                if jdata > 1
                    Ncells(2) = jdata -1;
                end
                Ncells(3) = kdata;
                if kdata > 1
                    Ncells(3) = kdata-1;
                end

                % load zone              
                if Prop.IVARLOC
                    ndata_cell = prod(Ncells)*sum(VARLOC==1);
                    ndata_nodal = ndata*sum(VARLOC==0);
                    ndata_total = ndata_nodal + ndata_cell;
                    % error('Don''t support IJK with varloc now');
                    if strcmpi(Prop.DATAPACKING,'POINT')
                        error('Don''t support IJK with varloc now with DATAPACKING=POINT')
                    end
                end

                    
                % let's load the data   
                
                if strcmpi(Prop.DATAPACKING,'POINT')
                    data_nodal = myfscanf(fid,'%f',ndata_total,IMY_FSCANF);
                    data_nodal = reshape(data_nodal,localNVAR,idata,jdata,kdata);
                    data2 = zeros(NVAR,idata,jdata,kdata);
                    if ~isempty(Prop.PASSIVEVARLIST) && ~nopassive
                        for i = 1:length(localVARindex)
                            data2(localVARindex(i),:,:,:) = data_nodal(i,:,:,:);
                        end
                    end
                elseif strcmpi(Prop.DATAPACKING,'BLOCK')
                    data = myfscanf(fid,'%f',ndata_total,IMY_FSCANF);

                    data_nodal = reshape(data(1:ndata_nodal),idata,jdata,kdata,sum(VARLOC==0));
                    data_cell = reshape(data(ndata_nodal+1:end),Ncells(1),Ncells(2),Ncells(3),sum(VARLOC==1));
                    % data = reshape(data,idata,jdata,kdata,localNVAR);
                    %% I'm going to ignore this for now
                    data2 = zeros(idata,jdata,kdata,NVAR);
                    if ~isempty(Prop.PASSIVEVARLIST) && ~nopassive
                        for i = 1:length(localVARindex)
                            data2(:,:,:,localVARindex(i)) = data(:,:,:,i);
                        end
                    end
                else
                    error('DATAPACKING not supported');
                end
                
                if ~isempty(Prop.PASSIVEVARLIST) && ~nopassive
                    %% I'm going to ignore this for now
                    data = squeeze(data2);
                else
                    % data = squeeze(data);
                    data_nodal = squeeze(data_nodal);
                    data_cell = squeeze(data_cell);
                end
                    
                if (jdata == 1 && kdata == 1)
                    % data = data'; %if x-y line data, transpose the data for easier reading
                    data_nodal = data_nodal'; %if x-y line data, transpose the data for easier reading
                    data_cell = data_cell'; %if x-y line data, transpose the data for easier reading
                end

                Propfields = fieldnames(Prop);
                for i =1:length(Propfields)
                    zone(nzone).(Propfields{i}) = Prop.(Propfields{i});
                end
                zone(nzone).data_cell  = data_cell;
                zone(nzone).data_nodal = data_nodal;
                if nopassive
                    zone(nzone).ZONEVARlist = ZONEVARlist;
                end
                
                clear Prop;

            elseif ((isfield(Prop,'N') && isfield(Prop,'E')) || ( isfield(Prop,'NODES') && isfield(Prop,'ELEMENTS')))
                % load finite element zone
                if (isfield(Prop,'N') && isfield(Prop,'E') )
                    % load zone information
                    Nnode = Prop.N;
                    Ecell = Prop.E;
                elseif ( isfield(Prop,'NODES') && isfield(Prop,'ELEMENTS'))
                    Nnode = Prop.NODES;
                    Ecell = Prop.ELEMENTS;
                else
                    error("No node information found");
                end
                
                if ~isfield(Prop,'DATAPACKING')
                    Prop.DATAPACKING = 'BLOCK';
                end
                if debug
                    fprintf('ZONE %d is a finite element zone\n',nzone);
                    fprintf('ZONE %d size: Node: %d Element: %d\n',nzone,Nnode,Ecell);
                end

                % start loading data
                data = repmat(struct('name',[],'data',[]),localNVAR,1);
                if strcmpi(Prop.DATAPACKING,'BLOCK')
                    for i = 1:localNVAR
                        data(i).name = ZONEVARlist{i};
                        ivarloc = Prop.VARLOC(localVARindex(i));
                        if ivarloc == 0
                            ndata = Nnode;
                        else
                            ndata = Ecell;
                        end
                        data(i).data = myfscanf(fid,'%e',ndata,IMY_FSCANF);
                    end
                elseif (strcmpi(Prop.DATAPACKING,'POINT') && ~Prop.IVARLOC)
                    %ensure not VARLOCATION and in point mode
                    data2 = myfscanf(fid,'%e',Nnode*localNVAR,IMY_FSCANF);
                    data2 = reshape(data2,localNVAR,Nnode);
                    
                    for i = 1:localNVAR
                        data(i).name = ZONEVARlist{i};
                        data(i).data = data2(i,:);
                    end
                end
                
                % loading connect information
                if strcmpi(Prop.ZONETYPE,'FEQUADRILATERAL')
                    conn = fscanf(fid,'%d',Ecell*4);
                    conn = reshape(conn,4,[]);
                elseif strcmpi(Prop.ZONETYPE,'FETRIANGLE')
                    conn = fscanf(fid,'%d',Ecell*3);
                    conn = reshape(conn,3,[]);
                elseif strcmpi(Prop.ZONETYPE,'FEBRICK')
                    conn = fscanf(fid,'%d',Ecell*8);
                    conn = reshape(conn,8,[]);
                elseif strcmpi(Prop.ZONETYPE,'FETETRAHEDRON')
                    conn = fscanf(fid,'%d',Ecell*4);
                    conn = reshape(conn,4,[]);   
                else
                    warning('FEPOLYGON and FEPOLYHEDRAL are not verfied')
                    conn = cell(1,Ecell);
                    for i = 1:Ecell
                        line = fgetl(fid);
                        conn{i} = sscanf(line,'%d');
                        conn{i} = conn{i}';
                    end
                end
                conn = conn';
                
                Propfields = fieldnames(Prop);
                for i =1:length(Propfields)
                    zone(nzone).(Propfields{i}) = Prop.(Propfields{i});
                end
                zone(nzone).data = data;
                zone(nzone).conn = conn;      
                zone(nzone).ZONEVARlist = ZONEVARlist;
                clear Prop;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ~isfield(Prop,'I')
                if ~isfield(Prop,'DATAPACKING')
                    Prop.DATAPACKING = 'POINT';
                end
                if strcmpi(Prop.DATAPACKING,'BLOCK')
                    error('assumed length data only supported for DATAPACKING=POINT')
                end
                
                % Assume still an IJK ordered ZONE
                if isfield(Prop,'ZONETYPE')
                    if ~strcmpi(Prop.ZONETYPE,'ORDERED')
                        error(['ZONE ',num2str(nzone),' type is not ordered']);
                    end
                else
                    Prop.ZONETYPE = 'ORDERED';
                end
                [data_nodal,N_data_tmp] = myfscanf(fid,repmat('%f',1,localNVAR),Inf,IMY_FSCANF);
                if mod(N_data_tmp,localNVAR)==0
                    % assume we can continue
                    Prop.I = N_data_tmp/localNVAR;
                    idata = Prop.I;
                    jdata = 1; kdata = 1;
                    if debug
                        fprintf('ZONE %d is a IJK ordered zone\n',nzone);
                        fprintf('ZONE %d Size: %d x %d x %d NVar: %d\n',nzone,idata,jdata,kdata,localNVAR);
                    end
                    ndata = idata*jdata*kdata;
                else
                    error('format not supported yet')
                end

                ndata_total = ndata*localNVAR;

                % load zone              
                if Prop.IVARLOC
                    error('Don''t support IJK with varloc now with DATAPACKING=POINT')
                end

                    
                % let's load the data
                data_nodal = reshape(data_nodal,localNVAR,idata,jdata,kdata);
                data_nodal = data_nodal'; %if x-y line data, transpose the data for easier reading

                Propfields = fieldnames(Prop);
                for i =1:length(Propfields)
                    zone(nzone).(Propfields{i}) = Prop.(Propfields{i});
                end
                zone(nzone).data_cell  = [];
                zone(nzone).data_nodal = data_nodal;
                if nopassive
                    zone(nzone).ZONEVARlist = ZONEVARlist;
                end
                
                clear Prop;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                error('format not supported yet')
            end
        end
    end
end            

zone = zone(1:nzone);
                                
fclose(fid);
if debug
    fprintf('===========================================\n');
end
end

function [list ] = strl2list(temp)
temp = strsplit(temp,',');
list = [];
for i = 1:length(temp)
    if ~contains(temp{i},'-')
        list = [list, str2num(temp{i})];
    else
        c = strrep(temp{i},'-',' ');
        c = str2num(c);
        list = [list, c(1):c(2)];
    end
end
end      

function [temp2,N] = myfscanf(fid,format,n,IMY_FSCANF)
%% A modified fscanf inorder to read fortran format scientific notation number
% like 1.234-100 = 1.234E-100
if IMY_FSCANF
    temp = textscan(fid,'%s',n,'CommentStyle','#');
    temp = temp{1};
    temp2 = regexprep(temp,'(-?)(\d+(.\d+)?|.\d+)([-+]\d+)','$1$2E$3');
    temp2 = cellfun(@str2double,temp2,'UniformOutput', false);
    temp2 = cell2mat(temp2);
    N = numel(temp2);
else
    [temp2,N] = fscanf(fid,format,n);
end
end