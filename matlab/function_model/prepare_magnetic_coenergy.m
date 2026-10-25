function data = prepare_magnetic_coenergy(dataset)
% Read saved research reconstructions, without fitting or running a study.
% dataset: 'Set2' (default), 'Set1', or 'pooled' (separate recording groups).
if nargin < 1
    dataset = 'Set2';
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
research = fullfile(root,'Resources','Research','magnetic_coenergy');
switch lower(char(dataset))
    case 'set1'
        snapshot = 'Matlab_200_201_202_FluxEval_snapshot.npz';
        results = fullfile('results_set1','reconstructions.npz');
        dataset = 'Set1';
    case 'set2'
        snapshot = 'data_snapshot.npz';
        results = fullfile('results','reconstructions.npz');
        dataset = 'Set2';
    case 'pooled'
        snapshot = 'pooled_snapshot.npz';
        results = fullfile('results_between_sets','pooled_reconstructions.npz');
        dataset = 'Set1 + Set2';
    otherwise
        error('Use Set1, Set2, or pooled.');
end
data.name = dataset;
data.source = fullfile(research,results);
data.snapshot = fullfile(research,snapshot);
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s'));
unzip(data.source,fullfile(folder,'results'));
unzip(data.snapshot,fullfile(folder,'snapshot'));
read_result = @(name) read_npy(fullfile(folder,'results',[name '.npy']));
read_snapshot = @(name) read_npy(fullfile(folder,'snapshot',[name '.npy']));
data.x = read_result('currents');
data.measured = read_result('measured_flux');
data.groups = read_result('groups')+1;
data.excitation = read_snapshot('excitation_ref');
data.speed = read_snapshot('speed');
data.weights = read_snapshot('weights');
metadata = jsondecode(read_snapshot('metadata'));
data.methods = {'Integration','RBF','Spline','Hybrid'};
for m = 1:numel(data.methods)
    method = data.methods{m};
    data.values{m} = [read_result([method '_W']),read_result([method '_psi'])];
    L = read_result([method '_L']);
    % Row-major matrix order: dd,dq,de,qd,qq,qe,ed,eq,ee.
    for j = 1:3
        for k = 1:3
            data.values{m}(:,end+1) = L(:,j,k);
        end
    end
end

% Do not show reference-dependent values when the reference curve is outside
% the measured 3D convex hull. The saved models otherwise allow continuation.
scale = max(data.x,[],1)-min(data.x,[],1);
hull = delaunayTriangulation(data.x./scale);
reference = [zeros(size(data.x,1),2),data.x(:,3)];
data.reference_supported = ~isnan(pointLocation(hull,reference./scale));
data.relative = [1 4 13];
for m = 1:numel(data.methods)
    data.values{m}(~data.reference_supported,data.relative) = nan;
end

ng = max(data.groups);
data.triangles = cell(ng,1);
data.labels = cell(ng,1);
data.mirror_residual = nan(size(data.measured));
data.mirror_ie_difference = nan(size(data.x,1),1);
for g = 1:ng
    rows = find(data.groups == g);
    x = data.x(rows,:);
    label = sprintf('%g A Soll | ie %.2f..%.2f A', ...
        data.excitation(rows(1)),min(x(:,3)),max(x(:,3)));
    if strcmp(dataset,'Set1 + Set2')
        source = metadata.groups(g).source;
        recording = regexp(source,'Set[12]','match','once');
        label = sprintf('%s | %s',recording,label);
    end
    data.labels{g} = label;
    % Preserve within-group gaps. Invalid values are removed from these
    % original triangles later, never reconnected across the missing region.
    [xy,first] = unique(x(:,1:2)./scale(1:2),'rows');
    if size(xy,1) < 3 || rank(xy-mean(xy,1)) < 2
        continue
    end
    dt = delaunayTriangulation(xy);
    tri = dt.ConnectivityList;
    edge = dt.edges;
    limit = 9*mean(sum((xy(edge(:,1),:)-xy(edge(:,2),:)).^2,2));
    keep = true(size(tri,1),1);
    for pair = [1 2;1 3;2 3]'
        keep = keep & sum((xy(tri(:,pair(1)),:)-xy(tri(:,pair(2)),:)).^2,2) < limit;
    end
    data.triangles{g} = rows(first(tri(keep,:)));
    % Existing q-parity hypotheses: psi_d even, psi_q odd. Compare rather
    % than enforce them, using only supported triangles of this recording.
    [cell_index,bary] = pointLocation(dt,[x(:,1),-x(:,2)]./scale(1:2));
    supported = isfinite(cell_index);
    supported(supported) = keep(cell_index(supported));
    vertices = rows(first(tri(cell_index(supported),:)));
    for c = 1:2
        field = data.measured(:,c);
        mirrored = sum(bary(supported,:).*field(vertices),2);
        data.mirror_residual(rows(supported),c) = ...
            field(rows(supported))-(-1)^(c-1)*mirrored;
    end
    ie = data.x(:,3);
    data.mirror_ie_difference(rows(supported)) = ...
        x(supported,3)-sum(bary(supported,:).*ie(vertices),2);
end
end

function value = read_npy(path)
% Numeric arrays are little-endian float64/int64; metadata is UTF-32 JSON.
% Reading these directly keeps this MATLAB entry point independent of Python.
fid = fopen(path,'r','ieee-le');
cleanup = onCleanup(@() fclose(fid));
fread(fid,6,'uint8');
version = fread(fid,2,'uint8');
if version(1) == 1
    count = fread(fid,1,'uint16');
else
    count = fread(fid,1,'uint32');
end
header = char(fread(fid,count,'*char')');
shape = regexp(header,'''shape'':\s*\(([^)]*)\)','tokens','once');
dims = sscanf(strrep(shape{1},',',' '),'%d')';
if contains(header,'''<f8''')
    value = fread(fid,prod(dims),'*double');
elseif contains(header,'''<i8''')
    value = fread(fid,prod(dims),'int64=>double');
elseif contains(header,'''<U')
    value = char(fread(fid,inf,'*uint32')');
    return
else
    error('Unsupported numeric NPY array: %s',path);
end
if numel(dims) > 1
    if contains(header,'''fortran_order'': True')
        value = reshape(value,dims);
    else
        value = permute(reshape(value,fliplr(dims)),numel(dims):-1:1);
    end
end
end
