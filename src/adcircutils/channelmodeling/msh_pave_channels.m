function m_new = msh_pave_channels(m, pave_opts)
pave_opts.out14 = sprintf('./meshes/ncv30/%s_%s', pave_opts.output_mesh_prefix, pave_opts.editname);
S = read_flowlinefile(pave_opts.flowlinefile, pave_opts.zshift, pave_opts.zshift_taper_elev_range, pave_opts.effective_width_ratio, pave_opts.max_width, pave_opts.min_width, false);
S = remove_duplicated_points_in_S(S,pave_opts.radius_to_merge_shppoints);
channel_spacing = ones(1,length(S))*pave_opts.resolution;
channel_spacing_end1 = channel_spacing;
channel_spacing_end2 = channel_spacing;
pave_opts.channel_spacing = channel_spacing;
pave_opts.channel_spacing_end1 = channel_spacing_end1;
pave_opts.channel_spacing_end2 = channel_spacing_end2;
m_new = pave_channel(m, pave_opts);
end

function area = triangle_area(x, y)
idx1 = [1 2 3];
idx2 = [2 3 1];
edgevec.x = x(:,idx1) - x(:,idx2);
edgevec.y = y(:,idx1) - y(:,idx2);
edgelen = sqrt(edgevec.x.^2 + edgevec.y.^2);
a = edgelen(:,1);
b = edgelen(:,2);
c = edgelen(:,3);
s = (a + b + c) / 2.0;
area = sqrt(s.*(s-a).*(s-b).*(s-c));
end

function min_dist = calculate_min_dist_seq_points(xy)
dxy = xy(1:end-1,:) - xy(2:end,:);
dists = dxy(:,1).^2 + dxy(:,2).^2;
min_dist = min(dists);
min_dist = sqrt(min_dist);
end

function distance = point_to_line_segment_distance(pxy, pxy1, pxy2)
    % pxy: coordinates of the point [px, py]
    % pxy1: coordinates of the first end of the line segment [px1, py1]
    % pxy2: coordinates of the second end of the line segment [px2, py2]

    % Calculate vector components of the line segment
    dx = pxy2(1) - pxy1(1);
    dy = pxy2(2) - pxy1(2);

    % Calculate squared distance from pxy1 to pxy2
    l_squared = dx^2 + dy^2;

    % Handle case where line segment is a single point
    if l_squared == 0
        distance = norm(pxy - pxy1);
        return;
    end

    % Calculate parameter t
    t = ((pxy(1) - pxy1(1)) * dx + (pxy(2) - pxy1(2)) * dy) / l_squared;

    % Clamp t to lie between 0 and 1
    t = max(0, min(1, t));

    % Calculate closest point on the line segment
    closest_point = [pxy1(1) + t * dx, pxy1(2) + t * dy];

    % Calculate distance between closest_point and pxy
    distance = norm(pxy - closest_point);
end

function esize = calculate_element_size(xy)
assert(all(size(xy) == [3, 2]))
dist1 = point_to_line_segment_distance(xy(1,:), xy(2,:), xy(3,:));
dist2 = point_to_line_segment_distance(xy(2,:), xy(3,:), xy(1,:));
dist3 = point_to_line_segment_distance(xy(3,:), xy(1,:), xy(2,:));
esize = min([dist1, dist2, dist3]);
end

function [ x,y ] = deg2cart( lat,lon,lat00,lon00 )
    %CPP_conv Converts lat and lon to x and y
    % lon0 =  mean(lon) * pi/180; lat0 = mean(lat) * pi/180;

    R = 6378206.4;
    lonr = lon * pi/180; latr = lat * pi/180;
    lon00r = lon00 * pi/180; lat00r = lat00 * pi/180;
    x = R * (lonr - lon00r) * cos(lat00r);
    y = R * latr;
end

function [ lat,lon ] = cart2deg( x,y,lat00,lon00 )
    %CPP_iconv Converts x and y to lat and lon
    
    R = 6378206.4;
    lon00r = lon00 * pi/180; lat00r = lat00 * pi/180;
    lon = rad2deg(lon00r + x/(R*cos(lat00r)));
    lat = rad2deg(y/R);
end

function compare_DEM_and_ADCIRC(dem,m)
dem_z_any = ones(size(m.p,1),1);
dem_z_any(:) = nan;
dem_z = ones(size(m.p,1),dem.ntiffs);
dem_z(:) = nan;
for i=1:size(m.p,1)
    lon = m.p(i,1);
    lat = m.p(i,2);
    for k=1:dem.ntiffs
        tifR = dem.tifRs{k};

        switch class(tifR)
            case 'map.rasterref.GeographicCellsReference'
                 [col, row] = geographicToIntrinsic(tifR, lat, lon);
        %          row = integer(row);
        %          col = integer(col);
                irow = round(row);
                icol = round(col);
                if irow < 1 || irow > size(dem.tifZs{k},1) || icol < 1 || icol > size(dem.tifZs{k},2)
                    continue
                end
                Vinterpolated = dem.tifZs{k}(irow,icol);
            case 'map.rasterref.MapCellsReference'
                [x,y] = projfwd(proj,lat,lon);
                if x < tifR.XWorldLimits(1) || x > tifR.XWorldLimits(2) || ...
                   y < tifR.YWorldLimits(1) || y > tifR.YWorldLimits(2)
                    continue
                end
                Vinterpolated = mapinterp(dem.tifZs{k},dem.tifRs{k},x,y);
            otherwise
                 error('Unrecognized format for R');
        end

        if ~isnan(Vinterpolated)
            dem_z_any(i) = -Vinterpolated;
            dem_z(i,k) = -Vinterpolated;
            break;
        end
    end
end
idx_any = ~isnan(dem_z_any);
idx = ~isnan(dem_z);
mean_dem = mean(-dem_z_any(idx_any));
mean_msh = mean(-m.b(idx_any));
fprintf('mean_dem: %.3f\n', mean_dem);
fprintf('mean_msh: %.3f\n', mean_msh);
fprintf('mean(dem-msh): %.3f\n', mean(-dem_z_any(idx_any) - -m.b(idx_any)));
fprintf('median(dem-msh): %.3f\n', median(-dem_z_any(idx_any) - -m.b(idx_any)));

figure
edges = -6:0.25:6;
for k=1:dem.ntiffs
    h = histogram(-dem_z(idx(:,k),k)- -m.b(idx(:,k)),edges);
    hold on
end

figure
plot([-100 100],[-100 100],'k-')
hold on
for k=1:dem.ntiffs
    plot(m.b(idx(:,k)),dem_z(idx(:,k),k),'.')
end
xlim([min([m.b(idx_any);dem_z_any(idx_any)]),max([m.b(idx_any);dem_z_any(idx_any)])])
ylim([min([m.b(idx_any);dem_z_any(idx_any)]),max([m.b(idx_any);dem_z_any(idx_any)])])
axis equal

m2 = m;
m2.b(:) = nan;
m2.b(idx_any) = -dem_z_any(idx_any) - -m.b(idx_any);
m2.write('dem_check_diff')

m2 = m;
m2.b(:) = nan;
m2.b(idx_any) = dem_z_any(idx_any);
m2.write('dem_check_demz')
end

function [mconn] = build_mesh_inside_polylines(polys_all, lat00, lon00, rseed, nnsub, msub, sub_nd_mesh_sizes, reject_margin, max_num_rejected, niter_relax, channel_spacing, nodes_inchannel, neinodes, ignore_depth_in_channel, no_channelmesh, plot_level)
% Generate nodes for the connecting mesh

[tmpx,tmpy] ...
    = deg2cart(msub.p(1,2),msub.p(1,1),lat00,lon00);

%write_edges(poly_inner_c,poly_outer_c);

offset = 0;
poly_nodes = [];
poly_edges = [];
for j=1:length(polys_all)
% for j=1
    poly = polys_all{j};
    poly = poly(1:end-1,:); % drop the last node assuming that it mathes the first node
    poly_nodes = [...
        poly_nodes; ...
        poly];
    ninnerj = size(poly,1);
    poly_edges = [...
        poly_edges; ...
        (offset+1:offset+ninnerj-1)' (offset+2:offset+ninnerj)'; ...
        offset+ninnerj offset+1];
    offset = offset + ninnerj;
end

if plot_level >= 2
    figure
    hold on
    for i=1:size(poly_edges,1)
        n1 = poly_edges(i,1);
        n2 = poly_edges(i,2);
        plot(poly_nodes([n1,n2],1),poly_nodes([n1,n2],2),'-');
        text(poly_nodes(n1,1), poly_nodes(n1,2), int2str(n1))
    end
    text(poly_nodes(n2,1), poly_nodes(n2,2), int2str(n2))
    axis equal
end

[new_nodes, new_elems] = mesh_gen(poly_nodes,poly_edges,channel_spacing,msub,sub_nd_mesh_sizes);

connect.x = new_nodes(:,1);
connect.y = new_nodes(:,2);

[connect.lat,connect.lon] = cart2deg(connect.x,connect.y,lat00,lon00);

% connect.lon(1:size(poly_nodes,1)) = [poly_outer(1:end-1,1);poly_inner_all(:,1)];
% connect.lat(1:size(poly_nodes,1)) = [poly_outer(1:end-1,2);poly_inner_all(:,2)];

connect.t = new_elems;

% Create the connecting mesh class
mconn = msh();
mconn.p = [connect.lon connect.lat];
mconn.t = connect.t;

% Give the connecting mesh the bathymetric depths
if ignore_depth_in_channel
    vals = make_vals_for_interp(nodes_inchannel,msub.b,neinodes);
else
    vals = msub.b;
end
mconn.b = interptri(msub.t,msub.p(:,1),msub.p(:,2),vals,mconn.p(:,1),mconn.p(:,2));

% Assign the depths at the nearest nodes if the above procedure did not
% assign valid values.
idx = find(isnan(mconn.b));
if ignore_depth_in_channel
    nodes_notinchannel = setdiff(1:size(msub.p,1),find(nodes_inchannel));
    ptmp = msub.p(nodes_notinchannel,:);
    jnd = dsearchn(ptmp,mconn.p(idx,:));
    mconn.b(idx) = msub.b(nodes_notinchannel(jnd));
else
    jnd = dsearchn(msub.p(:,:),double(msub.t),mconn.p(idx,:));
    mconn.b(idx) = msub.b(jnd);
end
end

function vals = make_vals_for_interp(nodes_inchannel,vals_org,neinodes)
vals = vals_org;
target_nodes = find(nodes_inchannel);
while ~isempty(target_nodes)
    target_nodes_next = [];
    for i = target_nodes'
        nei = neinodes{i};
        nei_not_inchannel = nei(~ismember(nei,find(nodes_inchannel)));
        if isempty(nei_not_inchannel)
            target_nodes_next(end+1) = i;
        else
            v = mean(vals_org(nei_not_inchannel));
            vals(i) = v;
        end
    end
    vals_org = vals;
    nodes_inchannel(:) = false;
    nodes_inchannel(target_nodes_next) = true;
    target_nodes = target_nodes_next';
end
end

function write_edges(poly_inner_c,poly_outer_c)
fid = fopen('edges.txt','w');
ninner = size(poly_inner_c,1)-1;
nouter = size(poly_outer_c,1)-1;
nnodes = ninner + nouter;
nedges = nnodes;

fprintf(fid,'%d %d\n',nnodes,nedges);

fprintf(fid,'%.3f %.3f\n',poly_inner_c(1:end-1,:)');
fprintf(fid,'%.3f %.3f\n',poly_outer_c(1:end-1,:)');

fprintf(fid,'%d %d\n',[(1:ninner-1)' (2:ninner)'; ninner 1]');
fprintf(fid,'%d %d\n',[(ninner+1:ninner+nouter-1)' (ninner+2:ninner+nouter)'; ninner+nouter ninner+1]');

fclose(fid);
end

function ma = add_condensed_nodes_to_f13(ma,condensed_nodes)
foundD = false;
for iAttrD = 1:length(ma.f13.defval.Atr)
    if strcmp('condensed_nodes',ma.f13.defval.Atr(iAttrD).AttrName)
        foundD = true;
        break;
    end
end

foundU = false;
for iAttrU = 1:length(ma.f13.userval.Atr)
    if strcmp('condensed_nodes',ma.f13.userval.Atr(iAttrU).AttrName)
        foundU = true;
        break;
    end
end

assert(foundD == foundU)

nval = size(condensed_nodes,2)-1;

if ~foundD
    ma = Calc_f13_Const(ma,'condensed_nodes',"1",nval,zeros(1,nval));
    iAttrD = length(ma.f13.userval.Atr);
    iAttrU = length(ma.f13.userval.Atr);
end

ValuesPerNode = ma.f13.defval.Atr(iAttrD).ValuesPerNode;
if ValuesPerNode < nval
    ma.f13.defval.Atr(iAttrD).ValuesPerNode = nval;
end
ma.f13.defval.Atr(iAttrD).Val(1:nval) = 0;

if  ~isempty(ma.f13.userval.Atr(iAttrU).Val)
    usernodes = ma.f13.userval.Atr(iAttrU).Val(1,:);
    assert(~any(ismember(usernodes,condensed_nodes(:,1)')))
end

ma.f13.userval.Atr(iAttrU).usernumnodes ...
    = ma.f13.userval.Atr(iAttrU).usernumnodes + size(condensed_nodes,1);

ma.f13.userval.Atr(iAttrU).Val(1:nval+1,end+1:end+size(condensed_nodes,1)) ...
    = condensed_nodes';
end

function ma = add_initial_river_elevation_to_f13(ma,channel_nodes)
foundD = false;
for iAttrD = 1:length(ma.f13.defval.Atr)
    if strcmp('initial_river_elevation',ma.f13.defval.Atr(iAttrD).AttrName)
        foundD = true;
        break;
    end
end

foundU = false;
for iAttrU = 1:length(ma.f13.userval.Atr)
    if strcmp('initial_river_elevation',ma.f13.userval.Atr(iAttrU).AttrName)
        foundU = true;
        break;
    end
end

assert(foundD == foundU)

nval = 1;
defval = -99999.0;

if ~foundD
    ma = Calc_f13_Const(ma,'initial_river_elevation',"m",nval,defval);
    iAttrD = length(ma.f13.userval.Atr);
    iAttrU = length(ma.f13.userval.Atr);
end

ValuesPerNode = ma.f13.defval.Atr(iAttrD).ValuesPerNode;
assert(ValuesPerNode == nval);
defval = ma.f13.defval.Atr(iAttrD).Val;

channelnodes = nonzeros(unique(channel_nodes(:)));

uservals_full = ones(1,size(ma.p,1))*defval;
if  ~isempty(ma.f13.userval.Atr(iAttrU).Val)
    usernodes = ma.f13.userval.Atr(iAttrU).Val(1,:);
    uservals = ma.f13.userval.Atr(iAttrU).Val(2,:);
    uservals_full(usernodes) = uservals;
end
uservals_full(channelnodes) = -ma.b(channelnodes) + 0.5; % 0.5 m above the bed elevation

defval_new = -99999.0;
usernodes_new = find(uservals_full ~= defval_new);
uservals_new = uservals_full(usernodes_new);

ma.f13.userval.Atr(iAttrU).usernumnodes ...
    = length(usernodes_new);

ma.f13.defval.Atr(iAttrD).Val = defval_new;

ma.f13.userval.Atr(iAttrU).Val ...
    = [usernodes_new; uservals_new];
end

function ma = set_channel_nodal_attributes(ma,channel_nodes,keydefvalues,keychannelvalues)
for attrname=keys(keydefvalues)'
    assert(isKey(keychannelvalues,attrname))
    defval = keydefvalues(attrname); defval = defval{1};
    channelval = keychannelvalues(attrname); channelval = channelval{1};

    foundD = false;
    for iAttrD = 1:length(ma.f13.defval.Atr)
        if strcmp(attrname,ma.f13.defval.Atr(iAttrD).AttrName)
            foundD = true;
            break;
        end
    end

    foundU = false;
    for iAttrU = 1:length(ma.f13.userval.Atr)
        if strcmp(attrname,ma.f13.userval.Atr(iAttrU).AttrName)
            foundU = true;
            break;
        end
    end
    
    assert(foundD == foundU)
    
    nval = length(defval);
    
    if ~foundD
        ma = Calc_f13_Const(ma,attrname,"m",nval,defval);
        iAttrD = length(ma.f13.defval.Atr);
        iAttrU = length(ma.f13.userval.Atr);
    end
    
    defval = ma.f13.defval.Atr(iAttrD).Val;
    ValuesPerNode = ma.f13.defval.Atr(iAttrD).ValuesPerNode;
    assert(ValuesPerNode == nval);
    
    channelnodes = nonzeros(unique(channel_nodes(:)));
    
    uservals_full = ones(nval,size(ma.p,1))*defval(1);
    if  ~isempty(ma.f13.userval.Atr(iAttrU).Val)
        usernodes = ma.f13.userval.Atr(iAttrU).Val(1,:);
        uservals = ma.f13.userval.Atr(iAttrU).Val(2:end,:);
        uservals_full(:,usernodes) = uservals;
    end
    uservals_full(:,channelnodes) = channelval(1);
    
    usernodes_new = find(uservals_full(1,:) ~= defval(1));
    uservals_new = uservals_full(usernodes_new);
    
    ma.f13.userval.Atr(iAttrU).usernumnodes ...
        = size(usernodes_new,2);
    
    ma.f13.userval.Atr(iAttrU).Val ...
        = [usernodes_new; uservals_new];
end
end

%
% Merge the close location points and set properties such as depth
%
function shp = merge_nearby_points(shp, radius_to_merge_shppoints)

% Initialize map
% pointmap = [];
% for i=1:shp.nseq
%     pointmap = [pointmap; shp.seq(i).seq' zeros(length(shp.seq(i).seq),1)];
% end
pointmap = (1:shp.npoint)';
pointmap(:,2) = 0;

% Merge
for i=1:shp.npoint
    % Find the points that have not been checked out
    nco.point = find(pointmap(:,2) == 0);
    nco.npoint = length(nco.point);
    if nco.npoint == 0
        break;
    end
    nco.lonlat = shp.lonlat(nco.point,:);
    
    % Compute the distance at all nodes from node i
    lon0 = shp.lonlat(i,1);
    lat0 = shp.lonlat(i,2);
    
    lons = [ones(1,nco.npoint)*lon0;nco.lonlat(:,1)'];
    lons = lons(:);
    lats = [ones(1,nco.npoint)*lat0;nco.lonlat(:,2)'];
    lats = lats(:);
    
    dists = m_lldist(lons,lats)*1000.0;
    dists = dists(1:2:end);
    
    % Find the same location nodes
    same_location_nodes = nco.point(dists < radius_to_merge_shppoints);
    
    % Register the same location nodes
    pointmap(same_location_nodes,2) = i;
end    

% Generate new point IDs by removing duplicated points
pt_remain = pointmap(:,1) == pointmap(:,2);
pm = pointmap(pointmap(:,1) == pointmap(:,2));
pointmap = [pointmap zeros(size(pointmap,1),1)];
pointmap(pointmap(:,1) == pointmap(:,2),3) = (1:size(pm,1))';
if length(shp.seq(1).depth) == 1
    % Make depth/width/spacing tables
    dps = zeros(1,shp.npoint);
    wds = zeros(1,shp.npoint);
    sps = zeros(1,shp.npoint);
    for i=1:shp.nseq
        dps(shp.seq(i).seq(1)) = shp.seq(i).depth;
        wds(shp.seq(i).seq(1)) = shp.seq(i).width;
        sps(shp.seq(i).seq(1)) = shp.seq(i).spacing;
        dps(shp.seq(i).seq(end)) = shp.seq(i).depth;
        wds(shp.seq(i).seq(end)) = shp.seq(i).width;
        sps(shp.seq(i).seq(end)) = shp.seq(i).spacing;
    end
    % Compute the mean depth/width/spacing at the merging points
    for i=1:shp.nseq
        shp.seq(i).depths = zeros(size(shp.seq(i).seq));
        shp.seq(i).widths = zeros(size(shp.seq(i).seq));
        shp.seq(i).spacings = zeros(size(shp.seq(i).seq));
        for j=[1 length(shp.seq(i).seq)]
            k = shp.seq(i).seq(j);
            mp = find(pointmap(:,2)==k);
            if length(mp) == 2  % take a mean value only at the merging points where two points are going to be merged.
                shp.seq(i).depths(j) = mean(dps(mp));
                shp.seq(i).widths(j) = mean(wds(mp));
                shp.seq(i).spacings(j) = mean(sps(mp));
            else
                shp.seq(i).depths(j) = shp.seq(i).depth;
                shp.seq(i).widths(j) = shp.seq(i).width;
                shp.seq(i).spacings(j) = shp.seq(i).spacing;
            end
        end
    end
    % Compute the weights of the mid points of each sequence
    for i=1:shp.nseq
        seq = shp.seq(i).seq;
        l = m_lldist(shp.lonlat(seq,1),shp.lonlat(seq,2))*1000.0;
        totl = sum(l);
        cuml = cumsum(l);
        shp.seq(i).weights = [0 cuml'/totl];
    end
    % Compute the depths/widths/spacings at the mid points
    for i=1:shp.nseq
        w = shp.seq(i).weights;
        dps = shp.seq(i).depths(1);
        dpe = shp.seq(i).depths(end);
        wds = shp.seq(i).widths(1);
        wde = shp.seq(i).widths(end);
        sps = shp.seq(i).spacings(1);
        spe = shp.seq(i).spacings(end);
        
        shp.seq(i).depths(2:end-1) = ...
            dps + w(2:end-1).*(dpe - dps);
        shp.seq(i).widths(2:end-1) = ...
            wds + w(2:end-1).*(wde - wds);
        shp.seq(i).spacings(2:end-1) = ...
            sps + w(2:end-1).*(spe - sps);
    end
end

% Create a new shp data after point merger
shp_org = shp;
clear shp;
shp.lonlat = shp_org.lonlat(pt_remain,:);
shp.coord = shp_org.coord(pt_remain,:);
if isfield(shp_org,'dists')
    shp.dists = shp_org.dists(pt_remain);
end
for i=1:shp_org.nseq
    shp.seq(i).seq = pointmap(pointmap(shp_org.seq(i).seq,2),3)';
end
if length(shp_org.seq(1).depth) == 1
    for i=1:shp_org.nseq
        shp.seq(i).depth = shp_org.seq(i).depths;
        shp.seq(i).width = shp_org.seq(i).widths;
        shp.seq(i).spacing = shp_org.seq(i).spacings;
    end
else
    for i=1:shp_org.nseq
        shp.seq(i).depth = shp_org.seq(i).depth;
        shp.seq(i).width = shp_org.seq(i).width;
        shp.seq(i).spacing = shp_org.seq(i).spacing;
    end
end
shp.npoint = size(shp.coord,1);
shp.nseq   = length(shp.seq);

end

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

function shp = clean_sequences(S, lat00, lon00, radius_to_merge_shppoints, diffuse_depth, plot_level)
%% Make point and node sequence arrays

% Make a table for all nodes
clear shp;
shp.lonlat = [];
shp.width = [];
shp.depth = [];
for i=1:length(S)
    jstart = size(shp.lonlat,1) + 1;
    shp.lonlat = [shp.lonlat; S(i).X(1:end) S(i).Y(1:end)];
    jend   = size(shp.lonlat,1);
    shp.seq(i).seq = jstart:jend;
    shp.seq(i).width = S(i).width;
    shp.seq(i).depth = S(i).depth;
    shp.seq(i).spacing = S(i).spacing;
end
[xs,ys] = deg2cart(shp.lonlat(:,2),shp.lonlat(:,1),lat00,lon00);
shp.coord = [xs ys];
shp.npoint = size(shp.lonlat,1);
shp.nseq   = length(shp.seq);

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(shp.seq)
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
    end
    disp(shp.nseq);
    axis equal
end
% Merge close location points
shp = merge_nearby_points(shp, radius_to_merge_shppoints);
% Plot for an intermediate check
% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(shp.seq)
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
        j1 = shp.seq(i).seq(1);
        j2 = shp.seq(i).seq(2);
        x = 0.5*(shp.coord(j1,1)+shp.coord(j2,1));
        y = 0.5*(shp.coord(j1,2)+shp.coord(j2,2));
        text(x,y,0,num2str(i))
    end
    axis equal
end


% Break up a sequence if an end point of another sequence is in the middle of the sequence

% Create an end point list
endpoints = [];
for i=1:shp.nseq
    endpoints = [endpoints shp.seq(i).seq(1) shp.seq(i).seq(end)];
end

% Break up
for i=1:shp.nseq
    seq = shp.seq(i).seq;
    breakpointflags = zeros(size(seq));
    for j=endpoints
        id = find(seq == j);
        breakpointflags(id) = 1;
        breakpointflags(1) = 0;
        breakpointflags(end) = 0;
    end
    breakpointids = find(breakpointflags==1);
    if ~isempty(breakpointids) > 0
        % Middle sub-sequence
        for k=1:length(breakpointids)-1
            iseq = length(shp.seq)+1;
            bp1 = breakpointids(k);
            bp2 = breakpointids(k+1);
            shp.seq(iseq).seq = shp.seq(i).seq(bp1:bp2);
            shp.seq(iseq).depth = shp.seq(i).depth(bp1:bp2);
            shp.seq(iseq).width = shp.seq(i).width(bp1:bp2);
            shp.seq(iseq).spacing = shp.seq(i).spacing(bp1:bp2);
        end
        % Last sub-sequence
        iseq = length(shp.seq)+1;
        bp1 = breakpointids(end);
        bp2 = length(seq);
        shp.seq(iseq).seq = shp.seq(i).seq(bp1:bp2);
        shp.seq(iseq).depth = shp.seq(i).depth(bp1:bp2);
        shp.seq(iseq).width = shp.seq(i).width(bp1:bp2);
        shp.seq(iseq).spacing = shp.seq(i).spacing(bp1:bp2);
        % First sub-sequence (overwriting an exisiting seq)
        subseq = seq(1:breakpointids(1));
        iseq = i;
        bp1 = 1;
        bp2 = breakpointids(1);
        shp.seq(iseq).seq = shp.seq(i).seq(bp1:bp2);
        shp.seq(iseq).depth = shp.seq(i).depth(bp1:bp2);
        shp.seq(iseq).width = shp.seq(i).width(bp1:bp2);
        shp.seq(iseq).spacing = shp.seq(i).spacing(bp1:bp2);
    end
end
shp.nseq = length(shp.seq);


% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(shp.seq)
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
        j1 = shp.seq(i).seq(1);
        j2 = shp.seq(i).seq(2);
        x = 0.5*(shp.coord(j1,1)+shp.coord(j2,1));
        y = 0.5*(shp.coord(j1,2)+shp.coord(j2,2));
        text(x,y,0,num2str(i))
    end
    for i=1:shp.npoint
        x = shp.coord(i,1);
        y = shp.coord(i,2);
        text(x,y,0,num2str(i),Color='g')
    end
    axis equal
end

% Merge pairs of sequences connected only with each other

% Create a table
endpoints = zeros(0,0);
for i=1:shp.nseq
    endpoints = [endpoints ...
        [shp.seq(i).seq(1) shp.seq(i).seq(end);i i]];
end

% Merge
endpoints_u = unique(endpoints(1,:));
newseqmap = 1:shp.nseq;
for i=endpoints_u
    dup = find(endpoints(1,:) == i);
    if length(dup) == 2
        iseqs = endpoints(2,dup);
        
        iseq = iseqs(1);
        while true
            if newseqmap(iseq) == iseq
                break;
            end
            iseq = newseqmap(iseq);
        end
        
        jseq = iseqs(2);
        while true
            if newseqmap(jseq) == jseq
                break;
            end
            jseq = newseqmap(jseq);
        end
        
        is = shp.seq(iseq).seq(1);
        ie = shp.seq(iseq).seq(end);
        js = shp.seq(jseq).seq(1);
        je = shp.seq(jseq).seq(end);
        
        seqleni = length(shp.seq(iseq).seq);
        seqlenj = length(shp.seq(jseq).seq);
        
        if is == js
            assert(is == i);
            inew1 = flip(1:seqleni);
            inew2 = 2:seqlenj;
            newseq = [shp.seq(iseq).seq(inew1) shp.seq(jseq).seq(inew2)];
            
            newdp = [shp.seq(iseq).depth(inew1) shp.seq(jseq).depth(inew2)];
            newwd = [shp.seq(iseq).width(inew1) shp.seq(jseq).width(inew2)];
            newsp = [shp.seq(iseq).spacing(inew1) shp.seq(jseq).spacing(inew2)];
        elseif is == je
            assert(is == i);
            inew1 = 1:seqlenj;
            inew2 = 2:seqleni;
            newseq = [shp.seq(jseq).seq(inew1) shp.seq(iseq).seq(inew2)];
            
            newdp = [shp.seq(jseq).depth(inew1) shp.seq(iseq).depth(inew2)];
            newwd = [shp.seq(jseq).width(inew1) shp.seq(iseq).width(inew2)];
            newsp = [shp.seq(jseq).spacing(inew1) shp.seq(iseq).spacing(inew2)];
        elseif ie == js
            assert(ie == i);
            inew1 = 1:seqleni;
            inew2 = 2:seqlenj;
            newseq = [shp.seq(iseq).seq(inew1) shp.seq(jseq).seq(inew2)];
            
            newdp = [shp.seq(iseq).depth(inew1) shp.seq(jseq).depth(inew2)];
            newwd = [shp.seq(iseq).width(inew1) shp.seq(jseq).width(inew2)];
            newsp = [shp.seq(iseq).spacing(inew1) shp.seq(jseq).spacing(inew2)];
        elseif ie == je
            assert(ie == i);
            inew1 = 1:seqleni;
            inew2 = flip(1:seqlenj-1);
            newseq = [shp.seq(iseq).seq(inew1) shp.seq(jseq).seq(inew2)];
            
            newdp = [shp.seq(iseq).depth(inew1) shp.seq(jseq).depth(inew2)];
            newwd = [shp.seq(iseq).width(inew1) shp.seq(jseq).width(inew2)];
            newsp = [shp.seq(iseq).spacing(inew1) shp.seq(jseq).spacing(inew2)];
        else
            assert(false);
        end
        shp.seq(iseq).seq = newseq;
        shp.seq(iseq).depth = newdp;
        shp.seq(iseq).width = newwd;
        shp.seq(iseq).spacing = newsp;
        
        shp.seq(jseq).seq = newseq;
        shp.seq(jseq).depth = newdp;
        shp.seq(jseq).width = newwd;
        shp.seq(jseq).spacing = newsp;
        
        newseqmap(jseq) = iseq;
    end
end
shp.seq(newseqmap ~= 1:shp.nseq) = [];
shp.nseq = length(shp.seq);

for i=1:shp.nseq
    if shp.seq(i).depth(end) == 0.39
        disp(i)
    end
end

% diffuse depth to nodes with nan depth
% Create a table
if diffuse_depth
    endpoints = zeros(0,0);
    for i=1:shp.nseq
        endpoints = [endpoints ...
            [shp.seq(i).seq(1) shp.seq(i).seq(end); i i]];
    end
    loopcnt = 0;
    anynan_cnt_prev = -1;
    anynan_cnt = 0;
    while true
        for i=1:shp.nseq
            if any(isnan(shp.seq(i).depth))
                seq = shp.seq(i);
                js = seq.seq(1);
                je = seq.seq(end);
                dups = find(endpoints(1,:) == js);
                dupe = find(endpoints(1,:) == je);
                assert(length(dups) >= 2 || length(dupe) >= 2)
                if length(dups) >= 2 && isnan(seq.depth(1))
                    iis = endpoints(2,dups);
                    for ii = iis
                        jj = find(shp.seq(ii).seq == js, 1);
                        seq.depth(1) = max(seq.depth(1), shp.seq(ii).depth(jj));
                    end
                end
                if length(dupe) >= 2 && isnan(seq.depth(end))
                    iis = endpoints(2,dupe);
                    for ii = iis
                        jj = find(shp.seq(ii).seq == je, 1);
                        seq.depth(end) = max(seq.depth(end), shp.seq(ii).depth(jj));
                    end
                end
                shp.seq(i).depth = seq.depth;
            end
        end
        for i=1:shp.nseq
            if any(isnan(shp.seq(i).depth))
                seq = shp.seq(i);
                js = seq.seq(1);
                je = seq.seq(end);
                dups = find(endpoints(1,:) == js);
                dupe = find(endpoints(1,:) == je);
                assert(length(dups) >= 2 || length(dupe) >= 2)
                if (length(dups) >= 2 && anynan_cnt ~= anynan_cnt_prev && isnan(seq.depth(1))) || all(isnan(seq.depth))
                    continue
                end
                if (length(dupe) >= 2 && anynan_cnt ~= anynan_cnt_prev && isnan(seq.depth(end))) || all(isnan(seq.depth))
                    continue
                end
                if (length(dups) == 1 || anynan_cnt == anynan_cnt_prev) && isnan(seq.depth(1))
                    jnotnan = find(~isnan(seq.depth),1,'first');
                    assert(~isempty(jnotnan))
                    seq.depth(1) = seq.depth(jnotnan);
                end
                if (length(dupe) == 1 || anynan_cnt == anynan_cnt_prev) && isnan(seq.depth(end))
                    jnotnan = find(~isnan(seq.depth),1,'last');
                    assert(~isempty(jnotnan))
                    seq.depth(end) = seq.depth(jnotnan);
                end
                jnan    = find(isnan(seq.depth));
                jnotnan = find(~isnan(seq.depth)) ;
                seq.depth(jnan) = interp1(double(jnotnan), seq.depth(jnotnan), double(jnan));
                shp.seq(i).depth = seq.depth;
            end
        end
        anynan_cnt_prev = anynan_cnt;
        anynan_cnt = 0;
        for i=1:shp.nseq
            % check
            if any(isnan(shp.seq(i).depth))
                anynan_cnt = anynan_cnt + 1;
            end
        end
        if anynan_cnt == 0
            break;
        end
        loopcnt = loopcnt + 1;
        if mod(loopcnt-1,100) == 0
            fprintf('diffusing depth ... loopcnt = %d, anynan_cnt = %d\n', loopcnt, anynan_cnt)
        end
        if loopcnt > 1000
            break
        end
    end
end

% Correction on main stream width
% - Create a table
endpoints = zeros(0,0);
for i=1:shp.nseq
    endpoints = [endpoints ...
        [shp.seq(i).seq(1) shp.seq(i).seq(end);i i]];
end

% - Correction on rapid expansion
for i = 1:shp.nseq
    if length(shp.seq(i).width) >= 3
        for j=1:length(shp.seq(i).width)-1
            shp.seq(i).width(j+1) = min(shp.seq(i).width(j+1),1.2*shp.seq(i).width(j));
        end
        for j=length(shp.seq(i).width):-1:2
            shp.seq(i).width(j-1) = min(shp.seq(i).width(j-1),1.2*shp.seq(i).width(j));
        end
    end
end

% - Correction for making the widths along the main stream the same.
idx = [2 3; 3 1; 1 2];
endpoints_u = unique(endpoints(1,:));
for i=endpoints_u
    dup = find(endpoints(1,:) == i);
    if length(dup) == 3
        iseqs = endpoints(2,dup);
        ws = [0.0 0.0 0.0];
        iposs = [0 0 0];
        for j=1:length(iseqs)
            iseq = iseqs(j);
            iposs(j) = find(shp.seq(iseq).seq == i);
            ws(j) = shp.seq(iseq).width(iposs(j));
        end
        [~, imin] = min(ws);
        iis = idx(imin,:);
        wmean = mean(ws(iis));
        shp.seq(iseqs(iis(1))).width(iposs(iis(1))) = wmean;
        shp.seq(iseqs(iis(2))).width(iposs(iis(2))) = wmean;
    end
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(shp.seq)
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
        j1 = shp.seq(i).seq(1);
        j2 = shp.seq(i).seq(2);
        x = 0.5*(shp.coord(j1,1)+shp.coord(j2,1));
        y = 0.5*(shp.coord(j1,2)+shp.coord(j2,2));
        text(x,y,0,num2str(i))
    end
    axis equal
end

end

function [nd_in, el_in, el_in_mem, msub, nnsub, subbars, subbarlens] ...
    = create_submesh_in_boundbox(S, m, channel_spacing, removal_range, plot_level)
% Bounding box of the channel mesh
xmin = min(S(1).X);
xmax = max(S(1).X);
ymin = min(S(1).Y);
ymax = max(S(1).Y);
if length(S) > 1
    for i=2:length(S)
        xmin = min([xmin,S(i).X]);
        xmax = max([xmax,S(i).X]);
        ymin = min([ymin,S(i).Y]);
        ymax = max([ymax,S(i).Y]);
    end
end
bboxch = [xmin xmax ymin ymax];
% Select elements near the node sequence by a box
% margin_deg_lon = (bboxch(2)-bboxch(1))*1.0; % 100% of width
% margin_deg_lat = (bboxch(4)-bboxch(3))*1.0; % 100% of lat
margin_deg_lon = channel_spacing*removal_range*2.5/10^5; % m to deg
margin_deg_lat = channel_spacing*removal_range*2.5/10^5; % m to deg
rect_in = [bboxch(1)-margin_deg_lon ...
    bboxch(2)+margin_deg_lon ...
    bboxch(3)-margin_deg_lat ...
    bboxch(4)+margin_deg_lat];
nd_in = find(m.p(:,1) >= rect_in(1) & m.p(:,1) <= rect_in(2) & m.p(:,2) >= rect_in(3) & m.p(:,2) <= rect_in(4));
el_in_mem1 = ismember(m.t(:,1),nd_in);
el_in_mem2 = ismember(m.t(:,2),nd_in);
el_in_mem3 = ismember(m.t(:,3),nd_in);
el_in_mem = el_in_mem1 | el_in_mem2 | el_in_mem3;
el_in = find(el_in_mem);
% Create a subset of the original mesh
msub = msh();
msub.t = m.t(el_in,:);
nodes_sub = sort(unique(msub.t(:)));
nd_in = nodes_sub;

map = [];
map(nodes_sub) = 1:length(nodes_sub);

msub.t = map(msub.t);
msub.p = m.p(nodes_sub,:);
msub.b = m.b(nodes_sub);
nnsub = size(msub.p,1);
nesub = size(msub.t,1);

% Remove nodal attributes at the removed nodes and remap
msub.f13 = m.f13;
for k=1:m.f13.nAttr
    nds = m.f13.userval.Atr(k).Val(1,:);
    ir = ismember(nds,nodes_sub);
    msub.f13.userval.Atr(k).Val(:,~ir) = [];
    msub.f13.userval.Atr(k).Val(1,:) = ...
        map(msub.f13.userval.Atr(k).Val(1,:));
    msub.f13.userval.Atr(k).usernumnodes ...
        = size(msub.f13.userval.Atr(k).Val,2);
end
% Get edge lengths
[subbars,subbarlens] = GetBarLengths(msub);
% Plot for an intermediate check
if plot_level >= 2
    figure;
    triplot(msub.t,msub.p(:,1),msub.p(:,2));
    axis equal
end
end

function [new_nodes,new_elems] = mesh_gen(nodes,edges,minh,msh,mshh)
%---------------------------------------------- do size-fun.
% olfs.dhdx = +0.08;
% 
% [vlfs,tlfs, ...
% hlfs] = lfshfn2(nodes,edges, ...
%                 []  ,olfs) ;
% 
% hlfs(hlfs<minh) = minh;

vlfs = msh.coord;
tlfs = msh.t;
hlfs = mshh';

[slfs] = idxtri2(vlfs,tlfs) ;

figure;
patch('faces',tlfs(:,1:3),'vertices',vlfs , ...
    'facevertexcdata' , hlfs, ...
    'facecolor','interp', ...
    'edgecolor','none') ;
hold on; axis image off;
clear title
title(['MESH-SIZE: KIND=DELAUNAY, |TRIA|=', ...
    num2str(size(tlfs,1))]) ;

%---------------------------------------------- do mesh-gen.
hfun = @trihfn2;
clear opts
% pave_opts.KIND = 'DELAUNAY';
opts.ref1 = 'preserve';
% opts.ref2 = 'preserve';
% opts.rho2 = 10;
% opts.siz1 = 10;
% opts.siz2 = 10;
opts.dtri = 'constrained';

[vert,etri,tria,tnum] ...
    = refine2_custom(nodes,edges,[],opts,hfun , ...
              vlfs,tlfs,slfs,hlfs);

figure;
patch('faces',tria(:,1:3),'vertices',vert, ...
    'facecolor','w', ...
    'edgecolor',[.2,.2,.2]) ;
hold on; axis image off;
patch('faces',edges(:,1:2),'vertices',nodes, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
title(['TRIA-MESH: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(tria,1))]) ;

drawnow;

set(figure(1),'units','normalized', ...
    'position',[.05,.50,.30,.35]) ;
set(figure(2),'units','normalized', ...
    'position',[.35,.50,.30,.35]) ;

new_nodes = vert;
new_elems = tria;

% %------------------------------------ remove nodes on edges
% nnodes = size(nnodes,1);
% nnewnodes = size(new_nodes,1);
% 
% noderemains = ones(nnewnodes,1);
% 
% for i = nnodes+1:nnewnodes
%     area = triangle_area()
% end

end

function [elem,exyrange] = find_element_containing_xy(coord,t,x,y,exyrange)
if isempty(exyrange)
    exmin = min([coord(t(:,1),1) coord(t(:,2),1) coord(t(:,3),1)],[],2);
    exmax = max([coord(t(:,1),1) coord(t(:,2),1) coord(t(:,3),1)],[],2);
    eymin = min([coord(t(:,1),2) coord(t(:,2),2) coord(t(:,3),2)],[],2);
    eymax = max([coord(t(:,1),2) coord(t(:,2),2) coord(t(:,3),2)],[],2);
    exyrange = [exmin,exmax,eymin,eymax];
else
    exmin = exyrange(:,1);
    exmax = exyrange(:,2);
    eymin = exyrange(:,3);
    eymax = exyrange(:,4);
end
target_elems = find(x>=exmin&x<=exmax&y>=eymin&y<=eymax);
% assert(~isempty(target_elems))
found = false;
for elem = target_elems'
    if is_inside_element(coord,t,elem,x,y)
        found = true;
        break;
    end
end
% assert(found)
if ~found
    elem = -1;
end
end

function inside = is_inside_element(coord,t,elem,x,y)
ns = t(elem,:);
n1 = ns(1); n2 = ns(2); n3 = ns(3);
tx0 = [coord(n1,1) coord(n2,1) coord(n3,1)];
tx1 = [x coord(n2,1) coord(n3,1)];
tx2 = [x coord(n3,1) coord(n1,1)];
tx3 = [x coord(n1,1) coord(n2,1)];
ty0 = [coord(n1,2) coord(n2,2) coord(n3,2)];
ty1 = [y coord(n2,2) coord(n3,2)];
ty2 = [y coord(n3,2) coord(n1,2)];
ty3 = [y coord(n1,2) coord(n2,2)];
ta0 = triangle_area(tx0,ty0);
ta1 = triangle_area(tx1,ty1);
ta2 = triangle_area(tx2,ty2);
ta3 = triangle_area(tx3,ty3);
if abs(ta0 - (ta1+ta2+ta3)) < 1e-5
    inside = true;
else
    inside = false;
end
end

function S = read_flowlinefile(flowlinefile, zshift, zshift_taper_elev_range, effective_width_ratio, max_width, min_width, plot_fig)

if strcmp(flowlinefile(end-7:end), '.geojson')
    S = read_geojsonfile(flowlinefile, zshift, zshift_taper_elev_range, effective_width_ratio, max_width, min_width);
else
    S = read_shpfile(flowlinefile, zshift, zshift_taper_elev_range, effective_width_ratio, max_width, min_width);
end

% replace -99999 with nan
for i = 1:length(S)
    % S(i).depth(S(i).depth==99999.0) = nan;
    % S(i).width(S(i).width==99999.0) = nan;
    S(i).depth(S(i).depth>20.0) = nan;
end

% Check duplication
duplicated = zeros(length(S),1);
for i = 1:length(S)
    for j = i+1:length(S)
        if length(S(i).X) ~= length(S(j).X)
            continue
        end
        dup = true;
        for k = 1:length(S(i).X)
            if S(i).X(k) ~= S(j).X(k) || S(i).Y(k) ~= S(j).Y(k)
                dup = false;
                break
            end
        end
        if dup
            if duplicated(i) == 0
                duplicated(i) = 1;
            end
            if duplicated(j) == 0
                duplicated(j) = -1;
            end
        end
    end
end

if any(duplicated)
    disp("duplicated lines found.")
    
    figure;
    hold on;
    for i=1:length(S)
        if duplicated(i)
            plot(S(i).X,S(i).Y,'ro-',LineWidth=10)
        else
            plot(S(i).X,S(i).Y,'k-')
        end
    end
    axis equal

    disp("removing one.")
    S(duplicated==-1) = [];
    % ME = MException("read_shpfile:inputError", "duplicated lines found.");
    % throw(ME)
    % return
end

% Plot for an intermediate check
if plot_fig
    figure;
    hold on;
    for i=1:length(S)
        plot(S(i).X,S(i).Y,'+-')
        text(0.5*(S(i).X(1)+S(i).X(2)),0.5*(S(i).Y(1)+S(i).Y(2)),0,num2str(i))
    end
    axis equal
end
end

function S = read_geojsonfile(geojsonfile, zshift, zshift_taper_elev_range, effective_width_ratio, max_width, min_width)
fid = fopen(geojsonfile);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
geoFlowline = jsondecode(str);

S = [];

nlines = length(geoFlowline.features);

for i=1:nlines
    S(i).X = geoFlowline.features(i).geometry.coordinates(:,1);
    S(i).Y = geoFlowline.features(i).geometry.coordinates(:,2);
    if size(S(i).X,1) ~= 1
        S(i).X = S(i).X';
        S(i).Y = S(i).Y';
    end

    if isfield(geoFlowline.features(i).properties,'depth')
        depth = geoFlowline.features(i).properties.depth;
        if isa(depth, 'char')
            depth = str2num(depth);
        end
        if ~isempty(depth)
            if isempty(zshift_taper_elev_range)
                zshift_rate = 1.0;
            else
                zshift_rate = 1 - (zshift_taper_elev_range(2) + depth)/(zshift_taper_elev_range(2) - zshift_taper_elev_range(1));
                zshift_rate = min(1.0,max(0.0,zshift_rate));
            end
            depth = depth + zshift_rate*zshift;
            S(i).depth = ones(size(S(i).X))*depth;
        end
    end
    if isfield(geoFlowline.features(i).properties,'pt_depth')
        pt_depth = geoFlowline.features(i).properties.pt_depth;
        if ~isempty(pt_depth)
            depths_str = split(pt_depth, ',');
            if length(depths_str) ~= length(S(i).X)
                disp([length(depths_str), length(S(i).X)])
                disp(depths_str)
                disp(S(i).X)
                disp('')
            end
            assert(length(depths_str) == length(S(i).X))
            depths = ones(1, length(depths_str));
            for k=1:length(depths)
                depth = double(string(depths_str(k))) * -1.0;
                if isempty(zshift_taper_elev_range)
                    zshift_rate = 1.0;
                else
                    zshift_rate = 1 - (zshift_taper_elev_range(2) + depth)/(zshift_taper_elev_range(2) - zshift_taper_elev_range(1));
                    zshift_rate = min(1.0,max(0.0,zshift_rate));
                end
                depths(k) = depth + zshift_rate * zshift;
            end
            S(i).depth = depths;
        end
    end
    found = false;
    if isfield(geoFlowline.features(i).properties,'width')
        width = geoFlowline.features(i).properties.width;
        found = true;
    end
    if isfield(geoFlowline.features(i).properties,'BtmWdth')
        width = geoFlowline.features(i).properties.BtmWdth;
        found = true;
    end
    if found
        if ~isempty(width)
            S(i).width = max(min_width, min(max_width, ones(size(S(i).X))*width*effective_width_ratio));
        end
    end
    if isfield(geoFlowline.features(i).properties,'pt_width')
        pt_width = geoFlowline.features(i).properties.pt_width;
        if ~isempty(pt_width)
            widths_str = split(pt_width, ',');
            assert(length(widths_str) == length(S(i).X))
            widths = ones(1, length(widths_str));
            for k=1:length(widths)
                widths(k) = double(string(widths_str(k)));
            end
            S(i).width = max(min_width, min(max_width, widths*effective_width_ratio));
        end
    end

    if ~isfield(S(i),'depth')
        S(i).depth = [];
    end
    if ~isfield(S(i),'width')
        S(i).width = [];
    end
end
end

function S = read_shpfile(shpfile, zshift, zshift_taper_elev_range, effective_width_ratio, max_width, min_width)
shpFlowline = m_shaperead(shpfile);
S = [];

for i=1:length(shpFlowline.ncst)
    S(i).X = shpFlowline.ncst{i}(:,1);
    S(i).Y = shpFlowline.ncst{i}(:,2);
    S(i).X(isnan(S(i).X)) = [];
    S(i).Y(isnan(S(i).Y)) = [];
    if size(S(i).X,1) ~= 1
        S(i).X = S(i).X';
        S(i).Y = S(i).Y';
    end

    for j=1:length(shpFlowline.fieldnames)
        if strcmp(shpFlowline.fieldnames{j},'depth')
            depth = shpFlowline.dbfdata{i,j};
            if ~isempty(depth)
                if isempty(zshift_taper_elev_range)
                    zshift_rate = 1.0;
                else
                    zshift_rate = 1 - (zshift_taper_elev_range(2) + depth)/(zshift_taper_elev_range(2) - zshift_taper_elev_range(1));
                    zshift_rate = min(1.0,max(0.0,zshift_rate));
                end
                depth = depth + zshift_rate * zshift;
                S(i).depth = ones(size(S(i).X))*depth;
            end
            break;
        end
    end
    for j=1:length(shpFlowline.fieldnames)
        if strcmp(shpFlowline.fieldnames{j},'pt_depth')
            pt_depth = shpFlowline.dbfdata{i,j};
            if ~isempty(pt_depth)
                depths_str = split(pt_depth, ',');
                if length(depths_str) ~= length(S(i).X)
                    disp([length(depths_str), length(S(i).X)])
                    disp(depths_str)
                    disp(S(i).X)
                    disp('')
                end
                assert(length(depths_str) == length(S(i).X))
                depths = ones(1, length(depths_str));
                for k=1:length(depths)
                    depth = double(string(depths_str(k)));
                    if isempty(zshift_taper_elev_range)
                        zshift_rate = 1.0;
                    else
                        zshift_rate = 1 - (zshift_taper_elev_range(2) + depth)/(zshift_taper_elev_range(2) - zshift_taper_elev_range(1));
                        zshift_rate = min(1.0,max(0.0,zshift_rate));
                    end
                    depths(k) = depth + zshift_rate * zshift;
                end
                S(i).depth = depths + zshift;
            end
            break;
        end
    end
    for j=1:length(shpFlowline.fieldnames)
        if strcmp(shpFlowline.fieldnames{j},'width') || ...
           strcmp(shpFlowline.fieldnames{j},'BtmWdth')
            width = shpFlowline.dbfdata{i,j};
            if ~isempty(width)
                S(i).width = max(min_width, min(max_width, ones(size(S(i).X))*width*effective_width_ratio));
            end
            break;
        end
    end
    for j=1:length(shpFlowline.fieldnames)
        if strcmp(shpFlowline.fieldnames{j},'pt_width')
            pt_width = shpFlowline.dbfdata{i,j};
            if ~isempty(pt_width)
                widths_str = split(pt_width, ',');
                assert(length(widths_str) == length(S(i).X))
                widths = ones(1, length(widths_str));
                for k=1:length(widths)
                    widths(k) = double(string(widths_str(k)));
                end
                S(i).width = max(min_width, min(max_width, widths*effective_width_ratio));
            end
            break;
        end
    end

    if ~isfield(S(i),'depth')
        S(i).depth = [];
    end
    if ~isfield(S(i),'width')
        S(i).width = [];
    end
end
end

function S = remove_duplicated_points_in_S(S, radius_to_merge_shppoints)
% Check duplicated points
for i=1:length(S)
    dists = m_lldist(S(i).X,S(i).Y)*1000.0;
    dup = dists < radius_to_merge_shppoints;
    S(i).X(dup) = [];
    S(i).Y(dup) = [];
end
end

function S = set_S_spacing(S, channel_spacing, channel_spacing_end1, channel_spacing_end2, lat00, lon00, link_with_width)
for i = 1:length(S)
    if link_with_width
        max_width = max(S(i).width);
    else
        max_width = 0.0;
    end
    S(i).spacing = max(max_width, channel_spacing(i));
    S(i).spacing_end1 = max(max_width, channel_spacing_end1(i));
    S(i).spacing_end2 = max(max_width, channel_spacing_end2(i));
end

figure
for i = 1:length(S)
    plot(S(i).X,S(i).Y,'+-')
    hold on
end
axis equal

end

function S = thinout_S(S, channel_spacing, channel_spacing_end1, channel_spacing_end2, lat00, lon00)
for i = 1:length(S)
for j = 1:length(S(i).X)
    lons = S(i).X;
    lats = S(i).Y;
    [xs,ys] = deg2cart(lats,lons,lat00,lon00);
    dl = channel_spacing(i);
    dl1 = channel_spacing_end1(i);
    dl2 = channel_spacing_end2(i);
    clnodex = xs(1);
    clnodey = ys(1);
    j = 1;
    i1 = 1;
    i2 = 2;
    while true
        x0 = clnodex(j);
        y0 = clnodey(j);
        x1 = xs(i1);
        y1 = ys(i1);
        x2 = xs(i2);
        y2 = ys(i2);
        dx1 = x1-x0;
        dy1 = y1-y0;
        dx2 = x2-x1;
        dy2 = y2-y1;
        a = dx2^2+dy2^2;
        b = dx1*dx2+dy1*dy2;
        c = dx1^2+dy1^2-dl^2;
        if(b^2 - a*c >= 0)
            tprime = (- b + sqrt(b^2 - a*c))/a;
        else
            tprime = 0.0;
        end

        if (tprime > 1.0) && i2 < length(xs)
            i1 = i1 + 1;
            i2 = i2 + 1;
            continue;
        end

        if (tprime < 0.0 || tprime > 1.0) && i2 == length(xs)
            break;
        end

        xp = x1 + tprime * (x2 - x1);
        yp = y1 + tprime * (y2 - y1);
        clnodex(end+1) = xp;
        clnodey(end+1) = yp;
        j = j + 1;        
    end
    clnodex(end+1) = xs(end);
    clnodey(end+1) = ys(end);

    [lats,lons] = cart2deg(clnodex,clnodey,lat00,lon00);
    S(i).X = lons;
    S(i).Y = lats;
    S(i).spacing = dl;
    S(i).spacing_end1 = dl1;
    S(i).spacing_end2 = dl2;
end
end

figure
for i = 1:length(S)
    plot(S(i).X,S(i).Y,'+-')
    hold on
end
axis equal

end

function shp = smooth_shp(shp, channel_smoothing_span)
    if channel_smoothing_span > 0
    for i=1:shp.nseq
%         l = length(shp.seq(i).width);
%         for j=1:l
%             j1 = max(1,j-5);
%             j2 = min(l,j+5);
%         end
%         shp.seq(i).width = smooth(shp.seq(i).width);
%         shp.seq(i).depth = smooth(shp.seq(i).depth);
        shp.seq(i).width = smoothdata(shp.seq(i).width,"rlowess",channel_smoothing_span);
        shp.seq(i).depth = smoothdata(shp.seq(i).depth,"rlowess",channel_smoothing_span);
    end
    end
end

function shp = thinout_shp(shp)
for iseq = 1:length(shp.seq)
    % Create a center line node sequence
    shpx = [cps.x;shp.coord(seq(2:end-1),1);cpe.x];
    shpy = [cps.y;shp.coord(seq(2:end-1),2);cpe.y];
    shpdp = shp.seq(k).depth;
    shpwd = shp.seq(k).width;
    shpdl= shp.dists(seq);
    clnodex = cps.x;
    clnodey = cps.y;
    clnodedp = cps.dp;
    clnodewd = cps.wd;
    clnodedl = shpdl(1);
    dl = shpdl(1);
    j = 1;
    i1 = 1;
    i2 = 2;
    while true
        x0 = clnodex(j);
        y0 = clnodey(j);
        dp0 = clnodedp(j);
        x1 = shpx(i1);
        y1 = shpy(i1);
        dp1 = shpdp(i1);
        wd1 = shpwd(i1);
        x2 = shpx(i2);
        y2 = shpy(i2);
        dp2 = shpdp(i2);
        wd2 = shpwd(i2);
        dx1 = x1-x0;
        dy1 = y1-y0;
        dx2 = x2-x1;
        dy2 = y2-y1;
        a = dx2^2+dy2^2;
        b = dx1*dx2+dy1*dy2;
        c = dx1^2+dy1^2-dl^2;
        if(b^2 - a*c >= 0)
            tprime = (- b + sqrt(b^2 - a*c))/a;
        else
            tprime = 0.0;
        end
        %         if (tprime < 0.0-1e-1 || tprime > 1.0+1e-1) && i2 < length(shpx)
        if (tprime > 1.0) && i2 < length(shpx)
            i1 = i1 + 1;
            i2 = i2 + 1;
            continue;
        end
        xp = x1 + tprime * (x2 - x1);
        yp = y1 + tprime * (y2 - y1);
        dpp = dp1 + tprime * (dp2 - dp1);
        wdp = wd1 + tprime * (wd2 - wd1);
        clnodex(end+1) = xp;
        clnodey(end+1) = yp;
        clnodedp(end+1) = dpp;
        clnodewd(end+1) = wdp;
        clnodedl(end+1) = dl;
        dl = shpdl(i1);
        j = j + 1;
        
        if (tprime < 0.0 || tprime > 1.0) && i2 == length(shpx)
            break;
        end
    end
    clnodex(end+1) = cpe.x;
    clnodey(end+1) = cpe.y;
    clnodedp(end+1) = cpe.dp;
    clnodewd(end+1) = cpe.wd;
    clnodedl(end+1) = dl;
    if length(clnodex) > 3
        x0 = 0.5*(clnodex(end-2) + clnodex(end));
        y0 = 0.5*(clnodey(end-2) + clnodey(end));
        dp0 = 0.5*(clnodedp(end-2) + clnodedp(end));
        wd0 = 0.5*(clnodewd(end-2) + clnodewd(end));
        ts = [];
        Ls = [];
        ii = [];
        for i=max(1,length(shpx)-2):length(shpx)-1
            ii(end+1) = i;
            x1 = shpx(i);
            y1 = shpy(i);
            dp1 = shpdp(i);
            wd1 = shpdp(i);
            x2 = shpx(i+1);
            y2 = shpy(i+1);
            dp2 = shpdp(i+1);
            Wd2 = shpwd(i+1);
            dx = x2 - x1;
            dy = y2 - y1;
            t = (dx*(x1-x0) + dy*(y1-y0))/(dx*dx + dy*dy);
            L = sqrt((x1+t*dx-x0)^2 + (y1+t*dy-y0)^2);
            ts(end+1) = t;
            Ls(end+1) = L;
        end
        Ls(ts < 0.0 | ts > 1.0) = 1e8;
        [M,I] = min(Ls);
        if ts(I) >= 0.0 && ts(I) <= 1.0 && Ls(I) ~= 1e8
            i = ii(I);
            x1 = shpx(i);
            y1 = shpy(i);
            dp1 = shpdp(i);
            wd1 = shpwd(i);
            x2 = shpx(i+1);
            y2 = shpy(i+1);
            dp2 = shpdp(i+1);
            wd2 = shpwd(i+1);
            dx = x2 - x1;
            dy = y2 - y1;
            ddp = dp2 - dp1;
            dwd = wd2 - wd1;
            t = ts(I);
            xt = x1 + t*dx;
            yt = y1 + t*dy;
            dpt = dp1 + t*ddp;
            wdt = wd1 + t*dwd;
        else
            xt = x0;
            yt = y0;
            dpt = dp0;
            wdt = wd0;
        end
        
        clnodex(end-1) = xt;
        clnodey(end-1) = yt;
        clnodedp(end-1) = dpt;
        clnodewd(end-1) = wdt;
    end
    while (sqrt((clnodex(end)-clnodex(end-1))^2+(clnodey(end)-clnodey(end-1))^2) ...
            < 0.5*clnodedl(end))
        if length(clnodex) == 2
            break
        end
        clnodex(end-1) = [];
        clnodey(end-1) = [];
        clnodedp(end-1) = [];
        clnodewd(end-1) = [];
        clnodedl(end-1) = [];
    end
    shp.seq(k).cl.coord = [clnodex' clnodey'];
    shp.seq(k).cl.depth = clnodedp;
    shp.seq(k).cl.width = clnodewd;
    shp.seq(k).cl.dists = clnodedl;
end
end

function plot_condensed_nodes(m)
grd.name = 'tmp';
grd.x = m.p(:,1);
grd.y = m.p(:,2);
grd.z = m.b;
grd.e = m.t;
grd.bnd = detbndy(grd.e);
nnm = size(m.p,1);

for i=1:m.f13.nAttr
    nvalpernd = m.f13.defval.Atr(i).ValuesPerNode;
    assert(nvalpernd >= 1);
    attrname_def = m.f13.defval.Atr(i).AttrName;
    attrname_usr = m.f13.userval.Atr(i).AttrName;
    assert(strcmp(attrname_def,attrname_usr));
    if ~strcmp(attrname_def,'condensed_nodes')
        continue;
    end
    vals = ones(nnm,1)*m.f13.defval.Atr(i).Val(1);
    idx = m.f13.userval.Atr(i).Val(:);
    idx(idx == 0) = [];
    vals(idx) = -1;
    vals(m.f13.userval.Atr(i).Val(1,:)) = m.f13.userval.Atr(i).Val(2,:);
    AXIS = [min(grd.x) max(grd.x) min(grd.y) max(grd.y)];
    CAXIS = [-1 1];
    if(CAXIS(1) == CAXIS(2))
        CAXIS = [CAXIS(1)-0.1 CAXIS(1)+0.1];
    end
    title = sprintf('%s', attrname_def);
    figure
    PPlotContour(grd,vals,AXIS,CAXIS,title);
end
end

function t = toc_disp(t,msg)
    disp(['Timer: ' msg])
    toc(t)
    t = tic;
end

function [tifZs, tifRs] = readCoNEDTiff(tifffiles,tiffdatumshifts,tiffdatumshiftlowerlimits)
tifZs = {};
tifRs = {};
for k=1:length(tifffiles)
    tifffile = tifffiles{k};
    fprintf('Reading tiff file %d/%d: %s\n', k, length(tifffiles), tifffile)
    [tifZs{k},tifRs{k}] = readgeoraster(tifffile);
    tifZs{k} = min(tifZs{k},max(tifZs{k} + tiffdatumshifts{k},tiffdatumshiftlowerlimits{k}));
end
end

%---------------------------------------------------------------
function m_updated_global = pave_channel(m, pave_opts)
%---------------------------------------------------------------
tStart = tic;
tPrev = tStart;

flowlinefile = pave_opts.flowlinefile;
dem = pave_opts.dem;
zshift = pave_opts.zshift;
zshift_taper_elev_range = pave_opts.zshift_taper_elev_range;
effective_width_ratio = pave_opts.effective_width_ratio;
max_width = pave_opts.max_width;
min_width = pave_opts.min_width;
channel_spacing = pave_opts.channel_spacing;
channel_spacing_end1 = pave_opts.channel_spacing_end1;
channel_spacing_end2 = pave_opts.channel_spacing_end2;
in14 = pave_opts.in14;
in13 = pave_opts.in13;
radius_to_merge_shppoints = pave_opts.radius_to_merge_shppoints;
removal_range = pave_opts.removal_range;
rep_dist_adjust = pave_opts.rep_dist_adjust;
rseed = pave_opts.rseed;
reject_margin = pave_opts.reject_margin;
max_num_rejected = pave_opts.max_num_rejected;
niter_relax = pave_opts.niter_relax;
ignore_depth_in_channel = pave_opts.ignore_depth_in_channel;
nudge_depth_at_channel_ends = pave_opts.nudge_depth_at_channel_ends;
no_channelmesh = pave_opts.no_channelmesh;
set_initial_river_elevation = pave_opts.set_initial_river_elevation;
default_ndattrs = pave_opts.default_ndattrs;
channel_ndattrs = pave_opts.channel_ndattrs;
channel_smoothing_span = pave_opts.channel_smoothing_span;
spacing_link_with_width = pave_opts.spacing_link_with_width;
min_spacing = pave_opts.min_spacing;
out14 = pave_opts.out14;
plot_level = pave_opts.plot_level;
write_output = pave_opts.write_output;
indomain_flowbc = pave_opts.indomain_flowbc;
do_renum = pave_opts.do_renum;

% width_default = 20.0;
% depth_default = 2.0;

S = read_flowlinefile( ...
    flowlinefile, zshift, zshift_taper_elev_range, ...
    effective_width_ratio, max_width, min_width, false);

lons = [];
lats = [];
for i=1:length(S)
    lons = [lons S(i).X];
    lats = [lats S(i).Y];
end
lon00 = mean(lons);
lat00 = mean(lats);

% S = thinout_S(S, channel_spacing, channel_spacing_end1, channel_spacing_end2, lat00, lon00);
S = set_S_spacing(S, channel_spacing, channel_spacing_end1, channel_spacing_end2, lat00, lon00, spacing_link_with_width);
S = remove_duplicated_points_in_S(S,radius_to_merge_shppoints);
% channel_spacing = ones(1,length(S))*channel_spacing; 
% channel_spacing_end1 = channel_spacing;
% channel_spacing_end2 = channel_spacing;

for i=1:length(S)
    S(i).len = length(S(i).X);
end

% % Checking if any string is isolated from the network.
% checked = zeros(length(S), 2);
% swept = zeros(length(S),1);
% i = 1; ii = 1;
% p = 1e-6;
% while true
%     checked(i,1) = 1;
%     checked(i,2) = 1;
%     swept(i) = 1;
%     found = false;
%     for iii = 1:S(i).len
%         lon0 = S(i).X(iii);
%         lat0 = S(i).Y(iii);
% 
%         for j=1:length(S)
%             if i == j
%                 continue
%             end
%             for jj = [1,2]
%                 if checked(j,jj) == 0
%                     if jj == 1
%                         lon1 = S(j).X(1);
%                         lat1 = S(j).Y(1);
%                     else
%                         lon1 = S(j).X(end);
%                         lat1 = S(j).Y(end);
%                     end
%                     if abs(lon0 - lon1) < p && abs(lat0 - lat1) < p
%                         checked(j,jj) = 1;
%                         found = true;
%                     end
%                 end
%             end
%         end
%     end
%     if ~found
%         for j=1:length(S)
%             if swept(j) == 0 && (checked(j,1) == 1 || checked(j,2) == 1)
%                 i = j;
%                 found = true;
%             end
%         end
%         if ~found
%             break
%         end
%     end
% end
% 
% % Plot for an intermediate check
% if plot_level >= 2 || any(swept ~= 1)
%     figure;
%     hold on;
% 
%     for i=1:length(S)
%         if swept(i) == 0
%             plot(S(i).X, S(i).Y, 'r-', LineWidth=10)
%         else
%             plot(S(i).X, S(i).Y, 'g-')
%         end
%     end
%     axis equal
% end
% 
% if any(swept ~= 1)
%     disp('There is at least one isolated string.')
%     return
% end

for i=1:length(S)
%     S(i).spacing = channel_spacing(i);    
    [x,y] = deg2cart(S(i).Y,S(i).X,lat00,lon00);
    S(i).X_cart = x;
    S(i).Y_cart = y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read and process finite element mesh 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the base mesh
if isempty(m)
    m = msh('fname',in14,'aux',{in13});
end
nn = size(m.p,1);
ne = size(m.t,1);

% Create a sub mesh around the sequances
[nd_in, el_in, el_in_mem, msub, nnsub, subbars, subbarlens] ...
    = create_submesh_in_boundbox(S, m, channel_spacing(1), removal_range, plot_level);
tPrev = toc_disp(tPrev,'1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve attributes from mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S2 = S;
% Genearate intermediate points in the sequences
for i=1:length(S)
    x = S(i).X_cart(1);
    y = S(i).Y_cart(1);
    cs = channel_spacing(i);
    S(i).spacing = ones(size(S(i).X_cart))*cs;
    s = S(i).spacing(1);
    if ~isempty(S(i).depth)
        d = S(i).depth(1);
    else
        d = [];
    end
    if ~isempty(S(i).width)
        w = S(i).width(1);
    else
        w = [];
    end
    for j=1:length(S(i).X_cart)-1
        x1 = S(i).X_cart(j);
        y1 = S(i).Y_cart(j);
        x2 = S(i).X_cart(j+1);
        y2 = S(i).Y_cart(j+1);
        s2 = S(i).spacing(j+1);
        if ~isempty(d)
            d2 = S(i).depth(j+1);
        end
        if ~isempty(w)
            w2 = S(i).width(j+1);
        end
        dx = x2 - x1;
        dy = y2 - y1;
        len = sqrt(dx*dx + dy*dy);
        num = floor(len / cs);
        if num <= 1
            x = [x;x2];
            y = [y;y2];
            s = [s;s2];
            if ~isempty(d)
                d = [d;d2];
            end
            if ~isempty(w)
                w = [w;w2];
            end
        else
            dlen = len / num;
            xs = x1 + dlen/len * dx * (1:num-1)';
            ys = y1 + dlen/len * dy * (1:num-1)';
            ss = ones(size(xs))*dlen;
            x = [x;xs;x2];
            y = [y;ys;y2];
            s = [s;ss;s2];
            if ~isempty(d)
                dd = ones(size(xs))*d2;
                d = [d;dd;d2];
            end
            if ~isempty(w)
                ww = ones(size(xs))*w2;
                w = [w;ww;w2];
            end
        end
    end
    [lat,lon] = cart2deg(x,y,lat00,lon00);
    S(i).X_cart = x;
    S(i).Y_cart = y;
    S(i).X = lon;
    S(i).Y = lat;
    S(i).spacing = s;
    if ~isempty(d)
        S(i).depth = d;
    end
    if ~isempty(w)
        S(i).width = w;
    end
    if size(S(i).spacing,1) ~= 1
        S(i).spacing = S(i).spacing';
    end
    if size(S(i).depth,1) ~= 1
        S(i).depth = S(i).depth';
    end
    if size(S(i).width,1) ~= 1
        S(i).width = S(i).width';
    end
end

% Apply channel_spacing_end1 and channel_spacing_end2
nendsecs = 10;
for i=1:length(S)
    cs_e1 = S(i).spacing_end1;
    cs_e2 = S(i).spacing_end2;
    npts = length(S(i).X_cart);
    for j=1:min(nendsecs,floor(npts/2))
        w = (j-1)/nendsecs;
        S(i).spacing(j) = (1-w)*cs_e1 + w*S(i).spacing(j);
        S(i).spacing(end-j+1) = (1-w)*cs_e2 + w*S(i).spacing(end-j+1);
    end
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(S)
        plot(S(i).X_cart,S(i).Y_cart,'+-')
    end
%     for i=1:length(S2)
%         plot(S2(i).X_cart,S2(i).Y_cart,'r+-')
%     end
    axis equal
end

% Compute the cartesian coordinates of msub nodes
[x,y] = deg2cart(msub.p(:,2),msub.p(:,1),lat00,lon00);
msub.coord = [x y];

tPrev = toc_disp(tPrev,'2');

% Make n2e
nnsub = size(msub.p,1);
nesub = size(msub.t,1);
n2esub = cell(nnsub,1);
for j = 1:nesub
    for nm = 1:3
        n2esub{msub.t(j,nm)} = [n2esub{msub.t(j,nm)} j];
    end
end
for i = 1:nnsub
    n2esub{i} = unique(n2esub{i});
end

% Make neinodes
neinodessub = cell(nnsub,1);
for i = 1:nnsub
    neielems = n2esub{i,1};
    nei = msub.t(neielems,:);
    neinodessub{i} = unique(nei(:)');
end

if write_output
    msub.write('msub');
end

tPrev = toc_disp(tPrev,'2');

msub_nodes_inchannel = [];
if no_channelmesh
    % Find elements and nodes along the shape lines
    msub_elems_inchannel = zeros(nesub,1);
    elems_at_ends = zeros(length(S),8);
    dist = 5.0;
    exyrange = [];
    for i=1:length(S)
        % tPrev = toc_disp(tPrev,'3.1');
        if mod(i-1,500) == 0
            fprintf('Finding elements along S(%d)\n',i)
        end
        for j=1:length(S(i).X_cart)
            x1 = S(i).X_cart(j);
            y1 = S(i).Y_cart(j);
        
            [elem,exyrange] = find_element_containing_xy(msub.coord,msub.t,x1,y1,exyrange);
    
            if elem > 0
                break;
            end
        end
        if elem <= 0 || j >= length(S(i).X_cart)
            continue;
        end
        assert(elem > 0)
        nm = msub.t(elem,:);
        
        msub_elems_inchannel(elem) = true;
    
        s_elems = elem;
    
        neielems = [];
        for neielem = n2esub(nm)'
            neielems = [neielems neielem{1}];
        end
        neielems = unique(neielems);
        
        elem_prev = elem;
    
        x0 = x1;
        y0 = y1;
        
        % tPrev = toc_disp(tPrev,'3.2');
    
        while true
            x2 = S(i).X_cart(j+1);
            y2 = S(i).Y_cart(j+1);
            
            dx = x2 - x0;
            dy = y2 - y0;
            dl = sqrt(dx*dx+dy*dy);
            
            dx1 = x1 - x0;
            dy1 = y1 - y0;
            dl1 = sqrt(dx1*dx1+dy1*dy1);
            
            if dl1 > dl
                j = j + 1;
                if j == length(S(i).X_cart)
                    break
                end
                x1 = S(i).X_cart(j);
                y1 = S(i).Y_cart(j);
                x0 = x1;
                y0 = y1;
                dist = 1.0;
                continue;
            end
            
            dx = dx/dl;
            dy = dy/dl;
            xp = x1 + dist * dx;
            yp = y1 + dist * dy;
            [elem,exyrange] = find_element_containing_xy(msub.coord,msub.t,xp,yp,exyrange);
            if elem <= 0
               break;
            end
            if elem == elem_prev
                x1 = xp;
                y1 = yp;
                %             fprintf('%.10f, %.10f\n',x1,y1)
                continue;
            else
                if ~any(ismember(elem,neielems))
                    dist = dist * 0.5;
                    if dist < 0.001
                        disp(dist)
                    end
                    assert(dist > 0.001);
                    continue;
                else
                    dist = dist * 2;
                end
                
                x1 = xp;
                y1 = yp;
                
                nm = msub.t(elem,:);
                
                neielems = [];
                for neielem = n2esub(nm)'
                    neielems = [neielems neielem{1}];
                end
                neielems = unique(neielems);
                
                msub_elems_inchannel(elem) = true;
    
                s_elems(end+1) = elem;
    
                elem_prev = elem;
            end
        end
    
        % tPrev = toc_disp(tPrev,'3.3');
    
        % store the elements found at the ends of this string
        if (length(s_elems) >= 8)
            elems_at_ends(i,1) = s_elems(1);
            if ~ismember(s_elems(2),elems_at_ends(i,:))
                elems_at_ends(i,2) = s_elems(2);
            end
            if ~ismember(s_elems(3),elems_at_ends(i,:))
                elems_at_ends(i,3) = s_elems(3);
            end
            if ~ismember(s_elems(4),elems_at_ends(i,:))
                elems_at_ends(i,4) = s_elems(4);
            end
            if ~ismember(s_elems(end-3),elems_at_ends(i,:))
                elems_at_ends(i,5) = s_elems(end-3);
            end
            if ~ismember(s_elems(end-2),elems_at_ends(i,:))
                elems_at_ends(i,6) = s_elems(end-2);
            end
            if ~ismember(s_elems(end-1),elems_at_ends(i,:))
                elems_at_ends(i,7) = s_elems(end-1);
            end
            if ~ismember(s_elems(end),elems_at_ends(i,:))
                elems_at_ends(i,8) = s_elems(end);
            end
        end
        % tPrev = toc_disp(tPrev,'3.4');
    
    end
    
    tPrev = toc_disp(tPrev,'3.5');
    
    % Remove the duplicated ends (elements at the branching points)
    elems_at_ends2 = elems_at_ends;
    for i = 1:length(S)
        elem_to_check1 = elems_at_ends(i,1);
        elem_to_check2 = elems_at_ends(i,2);
        elem_to_check3 = elems_at_ends(i,3);
        elem_to_check4 = elems_at_ends(i,4);
        S(i).is_end(1) = true;
        if elem_to_check1 ~= 0
            matched = find(elems_at_ends == elem_to_check1);
            if length(matched(:)) > 1
                elems_at_ends2(elems_at_ends2 == elem_to_check1) = 0;
                elems_at_ends2(elems_at_ends2 == elem_to_check2) = 0;
                elems_at_ends2(elems_at_ends2 == elem_to_check3) = 0;
                elems_at_ends2(elems_at_ends2 == elem_to_check4) = 0;
                S(i).is_end(1) = false;
            end
        end
        elem_to_check1 = elems_at_ends(i,8);
        elem_to_check2 = elems_at_ends(i,7);
        elem_to_check3 = elems_at_ends(i,6);
        elem_to_check4 = elems_at_ends(i,5);
        S(i).is_end(2) = true;
        if elem_to_check1 ~= 0
            matched = find(elems_at_ends == elem_to_check1);
            if length(matched(:)) > 1
                elems_at_ends2(elems_at_ends2 == elem_to_check1) = 0;
                elems_at_ends2(elems_at_ends2 == elem_to_check2) = 0;
                elems_at_ends2(elems_at_ends2 == elem_to_check3) = 0;
                elems_at_ends2(elems_at_ends2 == elem_to_check4) = 0;
                S(i).is_end(2) = false;
            end
        end
    end
    
    elems_at_ends = elems_at_ends2;
    
    tPrev = toc_disp(tPrev,'3');
    
    % Construct the list of the nodes in channels, not including ones at the ends 
    elems_at_ends = elems_at_ends(:);
    elems_at_ends(elems_at_ends == 0) = [];
    elems_inchannel = find(msub_elems_inchannel);
    if nudge_depth_at_channel_ends
        elems_inchannel = setdiff(elems_inchannel,elems_at_ends);
    end
    nodes_inchannel = msub.t(elems_inchannel,:);
    nodes_inchannel = unique(nodes_inchannel(:));
    msub_nodes_inchannel = zeros(nnsub,1);
    msub_nodes_inchannel(nodes_inchannel) = true;
    
    % Find elements containing each point
    % and set the width and depth to each shape points
    exyrange = [];
    for i=1:length(S)
        if isempty(S(i).width)
            S(i).width = ones(size(S(i).X_cart))*nan;
        end
        if isempty(S(i).depth)
            S(i).depth = ones(size(S(i).X_cart))*nan;
        end
        dp_interp = [];
        for j=1:length(S(i).X_cart)
            if ~isnan(S(i).width(j)) && ~isnan(S(i).depth(j))
                dp_interp(end+1) = S(i).depth(j);
                continue;
            end
            x = S(i).X_cart(j);
            y = S(i).Y_cart(j);
            [elem,exyrange] = find_element_containing_xy(msub.coord,msub.t,x,y,exyrange);
    %         assert(elem>0)    
            if elem <= 0
                % S(i).width(j) = nan;
                % dp_interp(end+1) = nan;
                % S(i).depth(j) = nan;
                dp_interp(end+1) = S(i).depth(j);
            else
                nids = [2 3;3 1;1 2];
                ns = msub.t(elem,:);
                xs = msub.coord(ns,1);
                ys = msub.coord(ns,2);
    
                if isnan(S(i).width(j))
                    nds = ns(nids);
                    dxs = msub.coord(nds(:,1),1) - msub.coord(nds(:,2),1);
                    dys = msub.coord(nds(:,1),2) - msub.coord(nds(:,2),2);
                    area = triangle_area(xs',ys');
                    dls = sqrt(dxs.*dxs+dys.*dys);
                    if j == 1
                        j1 = j;
                    else
                        j1 = j - 1;
                    end
                    if j == length(S(i).X_cart)
                        j2 = j;
                    else
                        j2 = j + 1;
                    end
                    assert(j1 ~= j2)
                    dx0 = S(i).X_cart(j2) - S(i).X_cart(j1);
                    dy0 = S(i).Y_cart(j2) - S(i).Y_cart(j1);
                    dl0 = sqrt(dx0*dx0+dy0*dy0);
                    dotps = dxs./dls*dx0/dl0 + dys./dls*dy0/dl0;
                    [vmax,imax] = max(abs(dotps));
                    len = sqrt(dxs(imax)*dxs(imax)+dys(imax)*dys(imax));
                    lenleg = area*2/len;
                    S(i).width(j) = min(lenleg,S(i).spacing(j));
                end
    
                if isnan(S(i).depth(j))
                    dp_interp(end+1) = interptri([1 2 3],msub.coord(ns,1),msub.coord(ns,2),msub.b(ns),x,y);
                    dp_max = max(msub.b(ns));
                    S(i).depth(j) = dp_max;
                else
                    dp_interp(end+1) = S(i).depth(j);
                end
            end
        end
        
        if S(i).is_end(1) && length(S(i).X_cart) >= 2
            S(i).depth(1) = dp_interp(1);
        end
        if S(i).is_end(2) && length(S(i).X_cart) >= 2
            S(i).depth(end) = dp_interp(end);
        end
    
        if size(S(i).width,1) ~= 1
            S(i).width = S(i).width';
        end
        if size(S(i).depth,1) ~= 1
            S(i).depth = S(i).depth';
        end
    end
end
% if ~no_channelmesh
    % Get depth values from DEM
    if isfield(dem,'ntiffs')
        for i=1:length(S)
            if mod(i-1,100) == 0
                fprintf('Processing DEM tiff for S(%d)\n', i)
            end
            for j=1:length(S(i).X_cart)
    %             if isnan(S(i).depth(j))
                    lon = S(i).X(j);
                    lat = S(i).Y(j);
                    for k=1:dem.ntiffs
                        tifR = dem.tifRs{k};
    
                        switch class(tifR)
                            case 'map.rasterref.GeographicCellsReference'
                                 [col, row] = geographicToIntrinsic(tifR, lat, lon);
                        %          row = integer(row);
                        %          col = integer(col);
                                irow = round(row);
                                icol = round(col);
                                if irow < 1 || irow > size(dem.tifZs{k},1) || icol < 1 || icol > size(dem.tifZs{k},2)
                                    continue
                                end
                                Vinterpolated = dem.tifZs{k}(irow,icol);
                            case 'map.rasterref.MapCellsReference'
                                proj = tifR.ProjectedCRS;
                                [x,y] = projfwd(proj,lat,lon);
                                if x < tifR.XWorldLimits(1) || x > tifR.XWorldLimits(2) || ...
                                   y < tifR.YWorldLimits(1) || y > tifR.YWorldLimits(2)
                                    continue
                                end
                                Vinterpolated = mapinterp(dem.tifZs{k},dem.tifRs{k},x,y);
                            otherwise
                                 error('Unrecognized format for R');
                        end
    
                        if ~isnan(Vinterpolated)
                            S(i).depth(j) = -Vinterpolated;
                            break;
                        end
                    end
    %             end
            end    
        end
    end
    % Replace NaN values with the non-NaN values found upstream or
    % downstream
    % for i=1:length(S)
    %     if ~any(isnan(S(i).width))
    %         continue
    %     end
    %     % - updatream (or downstream)
    %     dp_rep = nan;
    %     wd_rep = nan;
    %     for j=1:length(S(i).X_cart)
    %         if ~isnan(S(i).depth(j))
    %             dp_rep = S(i).depth(j);
    %         else
    %             S(i).depth(j) = dp_rep;
    %         end
    %         if ~isnan(S(i).width(j))
    %             wd_rep = S(i).width(j);
    %         else
    %             S(i).width(j) = wd_rep;
    %         end
    %     end
    %     % - downstream (or upstream)
    %     dp_rep = nan;
    %     wd_rep = nan;
    %     for j=flip(1:length(S(i).X_cart))
    %         if ~isnan(S(i).depth(j))
    %             dp_rep = S(i).depth(j);
    %         else
    %             S(i).depth(j) = dp_rep;
    %         end
    %         if ~isnan(S(i).width(j))
    %             wd_rep = S(i).width(j);
    %         else
    %             S(i).width(j) = wd_rep;
    %         end
    %     end
    % end
    % loopcnt = 0;
    % while true
    %     anynan_cnt = 0;
    %     for i=1:length(S)
    %         % depth
    %         if all(isnan(S(i).depth))
    %             for ii=1:length(S)
    %                 if ii==i
    %                     continue;
    %                 end
    %                 p = 1e-6;
    %                 if abs(S(i).X(1) - S(ii).X(1)) <= p && abs(S(i).Y(1) - S(ii).Y(1)) <= p
    %                     if ~isnan(S(ii).depth(1))
    %                         S(i).depth(:) = S(ii).depth(1);
    %                     end
    %                 elseif abs(S(i).X(1) - S(ii).X(end)) <= p && abs(S(i).Y(1) - S(ii).Y(end)) <= p
    %                     [i, ii] % 93 467
    %                     if ~isnan(S(ii).depth(end))
    %                         S(i).depth(:) = S(ii).depth(end);
    %                     end
    %                 elseif abs(S(i).X(end) - S(ii).X(1)) <= p && abs(S(i).Y(end) - S(ii).Y(1)) <= p
    %                     if ~isnan(S(ii).depth(1))
    %                         S(i).depth(:) = S(ii).depth(1);
    %                     end
    %                 elseif abs(S(i).X(end) - S(ii).X(end)) <= p && abs(S(i).Y(end) - S(ii).Y(end)) <= p
    %                     if ~isnan(S(ii).depth(end))
    %                         S(i).depth(:) = S(ii).depth(end);
    %                     end
    %                 end
    %             end
    %         end
    %         % width
    %         if all(isnan(S(i).width))
    %             for ii=1:length(S)
    %                 if ii==i
    %                     continue;
    %                 end
    %                 p = 1e-6;
    %                 if abs(S(i).X(1) - S(ii).X(1)) <= p && abs(S(i).Y(1) - S(ii).Y(1)) <= p
    %                     if ~isnan(S(ii).width(1))
    %                         S(i).width(:) = S(ii).width(1);
    %                     end
    %                 elseif abs(S(i).X(1) - S(ii).X(end)) <= p && abs(S(i).Y(1) - S(ii).Y(end)) <= p
    %                     if ~isnan(S(ii).width(end))
    %                         S(i).width(:) = S(ii).width(end);
    %                     end
    %                 elseif abs(S(i).X(end) - S(ii).X(1)) <= p && abs(S(i).Y(end) - S(ii).Y(1)) <= p
    %                     if ~isnan(S(ii).width(1))
    %                         S(i).width(:) = S(ii).width(1);
    %                     end
    %                 elseif abs(S(i).X(end) - S(ii).X(end)) <= p && abs(S(i).Y(end) - S(ii).Y(end)) <= p
    %                     if ~isnan(S(ii).width(end))
    %                         S(i).width(:) = S(ii).width(end);
    %                     end
    %                 end
    %             end
    %         end
    %         % check
    %         if any(isnan(S(i).width))
    %             anynan_cnt = anynan_cnt + 1;
    %         end
    %     end
    %     if anynan_cnt == 0
    %         break;
    %     end
    % 
    %     loopcnt = loopcnt + 1;
    %     if mod(loopcnt-1,100) == 0
    %         fprintf('diffusing depth/width ... loopcnt = %d\n', loopcnt)
    %     end
    % %     if loopcnt == 10
    % %         fprintf('Warning: forcing to set width and/or depth values.\n')
    % %         break;
    % %     end
    % end
% else
%     for i=1:length(S)
%         S(i).width = zeros(size(S(i).X));
%         S(i).depth = zeros(size(S(i).X));
%     end
% end
% % Force to set width and/or depth values.
% for i=1:length(S)
%     for j=1:length(S(i).width)
%         if isnan(S(i).width(j))
%             S(i).width(j) = 50.0;
%         end
%         if isnan(S(i).depth(j))
%             S(i).depth(j) = 0.0;
%         end
%     end
% end

% Plot for an intermediate check
if plot_level >= 2
    figure
    hold on
    for i=1:length(S)
        if any(isnan(S(i).width))
            plot(S(i).X,S(i).Y,'k-')
        else
            plot(S(i).X,S(i).Y,'-')
        end
    end
    axis equal
end

tPrev = toc_disp(tPrev,'4');

% Find nodes on the boundaries of msub mesh
boundary_nodes_sub = [];
boundary_elems_sub = [];
for i=1:nnsub
    neielems = n2esub{i};
    neinodes = neinodessub{i};
    neinodes(neinodes == i) = [];

    count = zeros(length(neinodes),1);
    for j=neielems
        nm = msub.t(j,:);
        i0 = find(nm == i);
        i1 = mod(i0-1+1,3)+1;
        i2 = mod(i0-1+2,3)+1;
        ii1 = find(neinodes == nm(i1));
        assert(~isempty(ii1));
        ii2 = find(neinodes == nm(i2));
        assert(~isempty(ii2));
        count(ii1) = count(ii1) + 1;
        count(ii2) = count(ii2) + 1;
    end

    if any(count == 1)
        boundary_nodes_sub(end+1) = i;
        boundary_elems_sub = ...
            [boundary_elems_sub neielems];
    end
end

boundary_nodes_sub = unique(boundary_nodes_sub);
% boundary_elems_sub = unique(boundary_elems_sub);
% nodes_tobekept_sub = msub.t(boundary_elems_sub,:);
nodes_tobekept_sub = boundary_nodes_sub;
nodes_tobekept_sub = unique(nodes_tobekept_sub(:));
nodes_tobekept = nd_in(nodes_tobekept_sub);

% Plot for an intermediate check
if plot_level >= 2
    figure;
    triplot(msub.t,msub.coord(:,1),msub.coord(:,2));
    hold on;
    plot(msub.coord(nodes_tobekept_sub,1),msub.coord(nodes_tobekept_sub,2),'o')
    plot(msub.coord(boundary_nodes_sub,1),msub.coord(boundary_nodes_sub,2),'o')
    axis equal
end

tPrev = toc_disp(tPrev,'5');

% Clean sequences
shp = clean_sequences(S, lat00, lon00, radius_to_merge_shppoints, pave_opts.diffuse_depth, plot_level);

plat = 36.467856;
plon = -77.635618;

knearest = 0;
dnearest = 99999.0;
for k = 1:shp.nseq
    seq = shp.seq(k);
    for i = 1:length(seq.seq)
        nd = seq.seq(i);
        ilat = shp.lonlat(nd,2);
        ilon = shp.lonlat(nd,1);
        d = sqrt((ilat-plat)^2 + (ilon-plon)^2);
        if d < dnearest
            knearest = k;
            dnearest = d;
        end
    end
end

ksnearest = 0;
dsnearest = 99999.0;
for k = 1:length(S)
    for i = 1:length(S(k).X)
        ilat = S(k).Y(i);
        ilon = S(k).X(i);
        d = sqrt((ilat-plat)^2 + (ilon-plon)^2);
        if d < dsnearest
            ksnearest = k;
            dsnearest = d;
        end
    end
end

figure
hold on
for k = 1:length(S)
    plot(S(k).X, S(k).Y, 'k-')
    if k == ksnearest
        plot(S(k).X, S(k).Y, 'r-')
    end
end
axis('equal')

figure
hold on
for k = 1:shp.nseq
    seq = shp.seq(k);
    plot(shp.lonlat(seq.seq,1), shp.lonlat(seq.seq,2), 'k-')
    if k == knearest
        plot(shp.lonlat(seq.seq,1), shp.lonlat(seq.seq,2), 'r-')
    end
end
axis('equal')

tPrev = toc_disp(tPrev,'6');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel mesh generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the representative mesh sizes at each shapefile node
nshpts = size(shp.lonlat,1);
shp.dists = zeros(nshpts,1);
count = zeros(nshpts,1);
for k=1:shp.nseq
    seq = shp.seq(k).seq;
    
    for j=1:length(seq)
        i = seq(j);
        shp.dists(i) = shp.dists(i) + shp.seq(k).spacing(j);
        count(i) = count(i) + 1;
    end
end

for i=1:nshpts
    shp.dists(i) = shp.dists(i) / count(i);
end

% Merge close location points
shp = merge_nearby_points(shp, radius_to_merge_shppoints);

% Apply smoothing
shp = smooth_shp(shp,channel_smoothing_span);

% Set spacing for short channels
for k=1:shp.nseq
    seq = shp.seq(k).seq;
    
    dxy = shp.coord(seq(1:end-1),:) - shp.coord(seq(2:end),:);
    seq_len = sum(sqrt(dxy(:,1).*dxy(:,1) + dxy(:,2).*dxy(:,2)));
    shp.seq(k).seq_len = seq_len;

    if length(seq) == 2
        shp.dists(seq(1)) = max(shp.dists(seq(1))*0.6, min(shp.dists(seq(1)), seq_len * 0.5));
        shp.dists(seq(2)) = max(shp.dists(seq(2))*0.6, min(shp.dists(seq(2)), seq_len * 0.5));
    end
end

% % Adjust spacing (thinning out)
% shp = thinout_shp(shp);

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(shp.seq)
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
    end
    axis equal
end

% Create an end point list
shp.endpoints.ids = [];
for i=1:shp.nseq
    assert(length(shp.seq(i).seq) >= 2);
    shp.endpoints.ids = [shp.endpoints.ids ...
        [shp.seq(i).seq(1) shp.seq(i).seq(end);
        shp.seq(i).seq(2) shp.seq(i).seq(end-1);
        i i;
        shp.seq(i).width(1) shp.seq(i).width(end);
        shp.seq(i).depth(1) shp.seq(i).depth(end);
        shp.seq(i).spacing(1) shp.seq(i).spacing(end)]];
end
shp.endpoints.uids = unique(shp.endpoints.ids(1,:));

% Find the connected points
for i=1:length(shp.endpoints.uids)
    i0 = shp.endpoints.uids(i);
    shp.endpoints.connected_points(i).ids = ...
        shp.endpoints.ids(2,shp.endpoints.ids(1,:) == i0);
    shp.endpoints.connected_points(i).seq = ...
        shp.endpoints.ids(3,shp.endpoints.ids(1,:) == i0);
    shp.endpoints.connected_points(i).width = ...
        shp.endpoints.ids(4,shp.endpoints.ids(1,:) == i0);
    shp.endpoints.connected_points(i).depth = ...
        shp.endpoints.ids(5,shp.endpoints.ids(1,:) == i0);
    shp.endpoints.connected_points(i).spacing = ...
        shp.endpoints.ids(6,shp.endpoints.ids(1,:) == i0);
    shp.endpoints.connected_points(i).depth(:) ...
        = max(shp.endpoints.connected_points(i).depth);
end

% Reorder
for i=1:length(shp.endpoints.uids)
    i0 = shp.endpoints.uids(i);
    shp.endpoints.connected_points(i).vecs = [];
    shp.endpoints.connected_points(i).lens = [];
    for j=1:length(shp.endpoints.connected_points(i).ids)
        i1 = shp.endpoints.connected_points(i).ids(j);
        vecx = shp.coord(i1,1) - shp.coord(i0,1);
        vecy = shp.coord(i1,2) - shp.coord(i0,2);
        len = sqrt(vecx*vecx + vecy*vecy);
        vecx = vecx/len;
        vecy = vecy/len;
        shp.endpoints.connected_points(i).vecs(end+1,1:3) = [vecx vecy 0];
        shp.endpoints.connected_points(i).lens(end+1) = len;
    end
    
    shp.endpoints.connected_points(i).rad = 0;
    vec1 = shp.endpoints.connected_points(i).vecs(1,:);
    for j=2:length(shp.endpoints.connected_points(i).ids)
        vec2 = shp.endpoints.connected_points(i).vecs(j,:);
        ac = acos(vec1*vec2');
        cp = cross(vec1,vec2);
        if cp(3) < 0
            ac = 2*pi-ac;
        end
        shp.endpoints.connected_points(i).rad(end+1) = ac;
    end
    
    if length(shp.endpoints.connected_points(i).ids) >= 3
        [sorted,I] = sort(shp.endpoints.connected_points(i).rad);
        shp.endpoints.connected_points(i).ids ...
            = shp.endpoints.connected_points(i).ids(I);
        shp.endpoints.connected_points(i).seq ...
            = shp.endpoints.connected_points(i).seq(I);
        shp.endpoints.connected_points(i).vecs ...
            = shp.endpoints.connected_points(i).vecs(I,:);
        shp.endpoints.connected_points(i).lens ...
            = shp.endpoints.connected_points(i).lens(I);
        shp.endpoints.connected_points(i).rad ...
            = shp.endpoints.connected_points(i).rad(I);
        shp.endpoints.connected_points(i).width ...
            = shp.endpoints.connected_points(i).width(I);
        shp.endpoints.connected_points(i).depth ...
            = shp.endpoints.connected_points(i).depth(I);
        shp.endpoints.connected_points(i).spacing ...
            = shp.endpoints.connected_points(i).spacing(I);
    end
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(shp.seq)
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
    end
    for i=1:length(shp.endpoints.uids)
        i0 = shp.endpoints.uids(i);
        for j=1:length(shp.endpoints.connected_points(i).ids)
            j1 = shp.endpoints.connected_points(i).ids(j);
            x = 0.5*(shp.coord(i0,1) + shp.coord(j1,1));
            y = 0.5*(shp.coord(i0,2) + shp.coord(j1,2));
            text(x,y,0,num2str(j))
        end
        i0 = shp.endpoints.uids(i);
        x0 = shp.coord(i0,1);
        y0 = shp.coord(i0,2);
        for j=1:length(shp.endpoints.connected_points(i).ids)
            x1 = x0 + shp.endpoints.connected_points(i).vecs(j,1)*1000.0;
            y1 = y0 + shp.endpoints.connected_points(i).vecs(j,2)*1000.0;
            plot([x0 x1],[y0 y1],'r-')
        end
    end
    axis equal
end

tPrev = toc_disp(tPrev,'7');

% Initialize a node id counter
gid_cnt = 1;
is_end_points = [];

% Determine the crossing points of neighboring branches
for i=1:length(shp.endpoints.uids)
    i0 = shp.endpoints.uids(i);
    x0 = shp.coord(i0,1);
    y0 = shp.coord(i0,2);
    
    assert(length(shp.endpoints.connected_points(i).ids) > 0);
    assert(length(shp.endpoints.connected_points(i).ids) ~= 2);
    
    if length(shp.endpoints.connected_points(i).ids) == 1
        vec1 = shp.endpoints.connected_points(i).vecs(1,:);
        w = shp.endpoints.connected_points(i).width(1);
        dp = shp.endpoints.connected_points(i).depth(1);
        sp = shp.endpoints.connected_points(i).spacing(1);
        nvec1 = 0.5*w*[-vec1(2) vec1(1)];
        p11x = x0 + nvec1(1);
        p11y = y0 + nvec1(2);
        p21x = x0 - nvec1(1);
        p21y = y0 - nvec1(2);
        shp.endpoints.connected_points(i).crossp(1).x = p11x;
        shp.endpoints.connected_points(i).crossp(1).y = p11y;
        shp.endpoints.connected_points(i).crossp(1).gid = gid_cnt;
        shp.endpoints.connected_points(i).crossp(1).depth = dp;
        shp.endpoints.connected_points(i).crossp(1).width = w;
        shp.endpoints.connected_points(i).crossp(1).spacing = sp;
        is_end_points(gid_cnt) = gid_cnt + 1;
        gid_cnt = gid_cnt + 1;
        shp.endpoints.connected_points(i).crossp(2).x = p21x;
        shp.endpoints.connected_points(i).crossp(2).y = p21y;
        shp.endpoints.connected_points(i).crossp(2).gid = gid_cnt;
        shp.endpoints.connected_points(i).crossp(2).depth = dp;
        shp.endpoints.connected_points(i).crossp(2).width = w;
        shp.endpoints.connected_points(i).crossp(2).spacing = sp;
        is_end_points(gid_cnt) =  -(gid_cnt - 1);
        gid_cnt = gid_cnt + 1;
    else
        for j=1:length(shp.endpoints.connected_points(i).ids)
            if j < length(shp.endpoints.connected_points(i).ids)
                jj = j+1;
            else
                jj = 1;
            end
            i1 = shp.endpoints.connected_points(i).ids(j);
            i2 = shp.endpoints.connected_points(i).ids(jj);
            vec1 = shp.endpoints.connected_points(i).vecs(j,:);
            vec2 = shp.endpoints.connected_points(i).vecs(jj,:);
            w1 = shp.endpoints.connected_points(i).width(j);
            w2 = shp.endpoints.connected_points(i).width(jj);
            dp1 = shp.endpoints.connected_points(i).depth(j);
            dp2 = shp.endpoints.connected_points(i).depth(jj);
            sp1 = shp.endpoints.connected_points(i).spacing(j);
            sp2 = shp.endpoints.connected_points(i).spacing(jj);
            if vec1(1)*vec2(1) + vec1(2)*vec2(2) < -0.90
                % fprintf('%d,%d,%d,%f,%f,%f\n',i0,i1,i2,vec1(1)*vec2(1) + vec1(2)*vec2(2),w1,w2)
                x1 = shp.coord(i1,1);
                y1 = shp.coord(i1,2);
                x2 = shp.coord(i2,1);
                y2 = shp.coord(i2,2);
                tx = x2 - x1;
                ty = y2 - y1;
                len = sqrt(tx*tx+ty*ty);
                nx = ty/len;
                ny = -tx/len;
                ln = 0.25*(w1+w2);
                p1x = x0 + nx*ln;
                p1y = y0 + ny*ln;
            else
                nvec1 = 0.5*w1*[-vec1(2) vec1(1)];
                nvec2 = 0.5*w2*[-vec2(2) vec2(1)];
                p11x = x0 + nvec1(1);
                p11y = y0 + nvec1(2);
                p21x = x0 - nvec2(1);
                p21y = y0 - nvec2(2);
                dett = -vec1(1)*vec2(2) + vec2(1)*vec1(2);
                px = -p11x + p21x;
                py = -p11y + p21y;
                t1 = (-vec2(2)*px + vec2(1)*py)/dett;
                p1x = p11x + t1*vec1(1);
                p1y = p11y + t1*vec1(2);
            end
            shp.endpoints.connected_points(i).crossp(j).x = p1x;
            shp.endpoints.connected_points(i).crossp(j).y = p1y;
            shp.endpoints.connected_points(i).crossp(j).gid = gid_cnt;
            shp.endpoints.connected_points(i).crossp(j).depth = 0.5*(dp1+dp2);
            shp.endpoints.connected_points(i).crossp(j).width = 0.5*(w1+w2);
            shp.endpoints.connected_points(i).crossp(j).spacing = 0.5*(sp1+sp2);
            is_end_points(gid_cnt) = 0;
            gid_cnt = gid_cnt + 1;
        end
    end
end
gid_endnode_max = gid_cnt-1;
assert(gid_endnode_max == length(is_end_points))

tPrev = toc_disp(tPrev,'8');

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(shp.seq)
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
    end
    for i=1:length(shp.endpoints.uids)
        if length(shp.endpoints.connected_points(i).ids) == 1
            nids = 2;
        else
            nids = length(shp.endpoints.connected_points(i).ids);
        end
        for j=1:nids
            x = shp.endpoints.connected_points(i).crossp(j).x;
            y = shp.endpoints.connected_points(i).crossp(j).y;
            gid = shp.endpoints.connected_points(i).crossp(j).gid;
            plot(x,y,'o')
            text(x,y,0,num2str(gid))
        end
    end
    axis equal
end

% Create center lines
for k=1:shp.nseq
    seq = shp.seq(k).seq;
    
    is = seq(1);
    ie = seq(end);
    
    is_ep = find(shp.endpoints.uids == is);
    ie_ep = find(shp.endpoints.uids == ie);
    
    is_cp = find(shp.endpoints.connected_points(is_ep).ids == seq(2));
    ie_cp = find(shp.endpoints.connected_points(ie_ep).ids == seq(end-1));
    
    is_crossp1 = is_cp;
    if is_cp == 1
        is_crossp2 = length(shp.endpoints.connected_points(is_ep).crossp);
    else
        is_crossp2 = is_cp-1;
    end
    
    ie_crossp1 = ie_cp;
    if ie_cp == 1
        ie_crossp2 = length(shp.endpoints.connected_points(ie_ep).crossp);
    else
        ie_crossp2 = ie_cp-1;
    end

    % Associate seq to crossp
    shp.seq(k).iconnected_points = [is_ep, ie_ep];
    shp.seq(k).icrossps = [is_crossp2, ie_crossp2];
    
    crossps1.x = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).x;
    crossps1.y = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).y;
    crossps1.dp = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).depth;
    crossps1.wd = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).width;
    crossps1.gid = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).gid;
    crossps2.x = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).x;
    crossps2.y = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).y;
    crossps2.dp = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).depth;
    crossps2.wd = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).width;
    crossps2.gid = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).gid;
    isendpoints = length(shp.endpoints.connected_points(is_ep).crossp) == 2;
    
    crosspe1.x = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).x;
    crosspe1.y = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).y;
    crosspe1.dp = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).depth;
    crosspe1.wd = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).width;
    crosspe1.gid = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).gid;
    crosspe2.x = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).x;
    crosspe2.y = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).y;
    crosspe2.dp = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).depth;
    crosspe2.wd = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).width;
    crosspe2.gid = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).gid;
    isendpointe = length(shp.endpoints.connected_points(ie_ep).crossp) == 2;
    
    nvecs = [crossps2.x-crossps1.x crossps2.y-crossps1.y];
    nvece = -[crosspe2.x-crosspe1.x crosspe2.y-crosspe1.y];
    lens = sqrt(nvecs*nvecs');
    lene = sqrt(nvece*nvece');
    nvecs = nvecs ./ lens;
    nvece = nvece ./ lene;
    tvecs = [-nvecs(2) nvecs(1)];
    tvece = [-nvece(2) nvece(1)];
        
    cps.x = 0.5*(crossps1.x + crossps2.x);
    cps.y = 0.5*(crossps1.y + crossps2.y);
    cps.dp = 0.5*(crossps1.dp + crossps2.dp);
    cps.wd = shp.seq(k).width(1);
    cpe.x = 0.5*(crosspe1.x + crosspe2.x);
    cpe.y = 0.5*(crosspe1.y + crosspe2.y);
    cpe.dp = 0.5*(crosspe1.dp + crosspe2.dp);
    cpe.wd = shp.seq(k).width(end);
    
    % Create a center line node sequence
    shpx = [cps.x;shp.coord(seq(2:end-1),1);cpe.x];
    shpy = [cps.y;shp.coord(seq(2:end-1),2);cpe.y];
    shpdp = shp.seq(k).depth;
    shpwd = shp.seq(k).width;
    shpdl= shp.dists(seq);
    clnodex = cps.x;
    clnodey = cps.y;
    clnodedp = cps.dp;
    clnodewd = cps.wd;
    clnodedl = shpdl(1);
    dl = shpdl(1)*1.1; % place the first node a little closer
    j = 1;
    i1 = 1;
    i2 = 2;
    while true
        x0 = clnodex(j);
        y0 = clnodey(j);
        dp0 = clnodedp(j);
        x1 = shpx(i1);
        y1 = shpy(i1);
        dp1 = shpdp(i1);
        wd1 = shpwd(i1);
        x2 = shpx(i2);
        y2 = shpy(i2);
        dp2 = shpdp(i2);
        wd2 = shpwd(i2);
        dx1 = x1-x0;
        dy1 = y1-y0;
        dx2 = x2-x1;
        dy2 = y2-y1;
        a = dx2^2+dy2^2;
        b = dx1*dx2+dy1*dy2;
        c = dx1^2+dy1^2-dl^2;
        if(b^2 - a*c >= 0)
            tprime = (- b + sqrt(b^2 - a*c))/a;
        else
            tprime = 0.0;
        end
        if (tprime > 1.0) && i2 < length(shpx)
            i1 = i1 + 1;
            i2 = i2 + 1;
            continue;
        end
        xp = x1 + tprime * (x2 - x1);
        yp = y1 + tprime * (y2 - y1);
        dpp = dp1 + tprime * (dp2 - dp1);
        wdp = wd1 + tprime * (wd2 - wd1);
        clnodex(end+1) = xp;
        clnodey(end+1) = yp;
        clnodedp(end+1) = dpp;
        clnodewd(end+1) = wdp;
        clnodedl(end+1) = dl;
        if i2 == length(shpx)
            dl = shpdl(i1) * 1.1; % 0.8;  % place the last nodes a little closer
        else
            dl = shpdl(i1);
        end
        j = j + 1;
        
        if i2 == length(shpx) && (tprime < 0.0 || tprime > 1.0)
            clnodex(end) = [];
            clnodey(end) = [];
            clnodedp(end) = [];
            clnodewd(end) = [];
            clnodedl(end) = [];
            break;
        end
    end

    clnodex(end+1) = cpe.x;
    clnodey(end+1) = cpe.y;
    clnodedp(end+1) = cpe.dp;
    clnodewd(end+1) = cpe.wd;
    clnodedl(end+1) = dl;

    if length(clnodex) >= 3
        if (clnodex(end-1) - clnodex(end))*(clnodex(end-2) - clnodex(end)) + (clnodey(end-1) - clnodey(end))*(clnodey(end-2) - clnodey(end)) <= 0.0
            clnodex(end-1) = [];
            clnodey(end-1) = [];
            clnodedp(end-1) = [];
            clnodewd(end-1) = [];
            clnodedl(end-1) = [];
        end
    end

    while true
        if length(clnodex) >= 3
            dxy = [clnodex(end) - clnodex(end-1), clnodey(end) - clnodey(end-1)];
            if norm(dxy) <= min_spacing
                clnodex(end-1) = [];
                clnodey(end-1) = [];
                clnodedp(end-1) = [];
                clnodewd(end-1) = [];
                clnodedl(end-1) = [];
            else
                break
            end
        else
            break
        end
    end

    % if length(clnodex) > 3
    %     x0 = 0.5*(clnodex(end-2) + clnodex(end));
    %     y0 = 0.5*(clnodey(end-2) + clnodey(end));
    %     dp0 = 0.5*(clnodedp(end-2) + clnodedp(end));
    %     wd0 = 0.5*(clnodewd(end-2) + clnodewd(end));
    %     ts = [];
    %     Ls = [];
    %     ii = [];
    %     for i=max(1,length(shpx)-2):length(shpx)-1
    %         ii(end+1) = i;
    %         x1 = shpx(i);
    %         y1 = shpy(i);
    %         dp1 = shpdp(i);
    %         wd1 = shpdp(i);
    %         x2 = shpx(i+1);
    %         y2 = shpy(i+1);
    %         dp2 = shpdp(i+1);
    %         Wd2 = shpwd(i+1);
    %         dx = x2 - x1;
    %         dy = y2 - y1;
    %         t = (dx*(x1-x0) + dy*(y1-y0))/(dx*dx + dy*dy);
    %         L = sqrt((x1+t*dx-x0)^2 + (y1+t*dy-y0)^2);
    %         ts(end+1) = t;
    %         Ls(end+1) = L;
    %     end
    %     Ls(ts < 0.0 | ts > 1.0) = 1e8;
    %     [M,I] = min(Ls);
    %     if ts(I) >= 0.0 && ts(I) <= 1.0 && Ls(I) ~= 1e8
    %         i = ii(I);
    %         x1 = shpx(i);
    %         y1 = shpy(i);
    %         dp1 = shpdp(i);
    %         wd1 = shpwd(i);
    %         x2 = shpx(i+1);
    %         y2 = shpy(i+1);
    %         dp2 = shpdp(i+1);
    %         wd2 = shpwd(i+1);
    %         dx = x2 - x1;
    %         dy = y2 - y1;
    %         ddp = dp2 - dp1;
    %         dwd = wd2 - wd1;
    %         t = ts(I);
    %         xt = x1 + t*dx;
    %         yt = y1 + t*dy;
    %         dpt = dp1 + t*ddp;
    %         wdt = wd1 + t*dwd;
    %     else
    %         xt = x0;
    %         yt = y0;
    %         dpt = dp0;
    %         wdt = wd0;
    %     end
    % 
    %     clnodex(end-1) = xt;
    %     clnodey(end-1) = yt;
    %     clnodedp(end-1) = dpt;
    %     clnodewd(end-1) = wdt;
    % end
    % while (sqrt((clnodex(end)-clnodex(end-1))^2+(clnodey(end)-clnodey(end-1))^2) ...
    %         < 0.5*clnodedl(end))
    %     if length(clnodex) == 2
    %         break
    %     end
    %     clnodex(end-1) = [];
    %     clnodey(end-1) = [];
    %     clnodedp(end-1) = [];
    %     clnodewd(end-1) = [];
    %     clnodedl(end-1) = [];
    % end
    shp.seq(k).cl.coord = [clnodex' clnodey'];
    shp.seq(k).cl.depth = clnodedp;
    shp.seq(k).cl.width = clnodewd;
    shp.seq(k).cl.dists = clnodedl;

    % Update direction vectors at endpoints
    dx = clnodex(2) - clnodex(1);
    dy = clnodey(2) - clnodey(1);
    len = sqrt(dx*dx + dy*dy);
    shp.endpoints.connected_points(is_ep).vecs(is_cp,:) = [dx/len, dy/len, 0];

    dx = clnodex(end-1) - clnodex(end);
    dy = clnodey(end-1) - clnodey(end);
    len = sqrt(dx*dx + dy*dy);
    shp.endpoints.connected_points(ie_ep).vecs(ie_cp,:) = [dx/len, dy/len, 0];
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(S)
        plot(S(i).X_cart,S(i).Y_cart,'k.-')
    end
    for i=1:shp.nseq
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
        plot(shp.seq(i).cl.coord(:,1),shp.seq(i).cl.coord(:,2),'*')
        for j=1:size(shp.seq(i).cl.coord,1)
            text(shp.seq(i).cl.coord(j,1), shp.seq(i).cl.coord(j,2), int2str(j))
        end
    end
    axis equal
end

% Associate seq.cl to endpoints
for i=1:shp.nseq
    shp.seq(i).cl.endpoints = [0 0];
end
for i=1:length(shp.endpoints.connected_points)
    for j=1:length(shp.endpoints.connected_points(i).seq)
        iseq = shp.endpoints.connected_points(i).seq(j);
        id = shp.endpoints.connected_points(i).ids(j);
        if shp.seq(iseq).seq(2) == id
            shp.seq(iseq).cl.endpoints(1) = i;
        elseif shp.seq(iseq).seq(end-1) == id
            shp.seq(iseq).cl.endpoints(2) = i;
        else
            assert(0);
        end
    end
end


if false
    figure
    hold on;
    for i=1:shp.nseq
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
    end
end

% - Update the crossing points of neighboring branches
for i=1:length(shp.endpoints.uids)
    i0 = shp.endpoints.uids(i);
    x0 = shp.coord(i0,1);
    y0 = shp.coord(i0,2);
    
    assert(length(shp.endpoints.connected_points(i).ids) > 0);
    assert(length(shp.endpoints.connected_points(i).ids) ~= 2);
    
    if length(shp.endpoints.connected_points(i).ids) == 1
        % Do nothing here
    else
        for j=1:length(shp.endpoints.connected_points(i).ids)
            if j < length(shp.endpoints.connected_points(i).ids)
                jj = j+1;
            else
                jj = 1;
            end
            i1 = shp.endpoints.connected_points(i).ids(j);
            i2 = shp.endpoints.connected_points(i).ids(jj);
            vec1 = shp.endpoints.connected_points(i).vecs(j,:);
            vec2 = shp.endpoints.connected_points(i).vecs(jj,:);
            w1 = shp.endpoints.connected_points(i).width(j);
            w2 = shp.endpoints.connected_points(i).width(jj);
            if vec1(1)*vec2(1) + vec1(2)*vec2(2) < -0.90
                % fprintf('%d,%d,%d,%f,%f,%f\n',i0,i1,i2,vec1(1)*vec2(1) + vec1(2)*vec2(2),w1,w2)
                x1 = shp.coord(i1,1);
                y1 = shp.coord(i1,2);
                x2 = shp.coord(i2,1);
                y2 = shp.coord(i2,2);
                tx = x2 - x1;
                ty = y2 - y1;
                len = sqrt(tx*tx+ty*ty);
                nx = ty/len;
                ny = -tx/len;
                ln = 0.25*(w1+w2);
                p1x = x0 + nx*ln;
                p1y = y0 + ny*ln;
            else
                nvec1 = 0.5*w1*[-vec1(2) vec1(1)];
                nvec2 = 0.5*w2*[-vec2(2) vec2(1)];
                p11x = x0 + nvec1(1);
                p11y = y0 + nvec1(2);
                p21x = x0 - nvec2(1);
                p21y = y0 - nvec2(2);
                dett = -vec1(1)*vec2(2) + vec2(1)*vec1(2);
                px = -p11x + p21x;
                py = -p11y + p21y;
                t1 = (-vec2(2)*px + vec2(1)*py)/dett;
                p1x = p11x + t1*vec1(1);
                p1y = p11y + t1*vec1(2);
            end
            shp.endpoints.connected_points(i).crossp(j).x = p1x;
            shp.endpoints.connected_points(i).crossp(j).y = p1y;
        end
    end
end

% - Create bank nodes
for k=1:shp.nseq
    seq = shp.seq(k).seq;
    
    is = seq(1);
    ie = seq(end);
    
    is_ep = find(shp.endpoints.uids == is);
    ie_ep = find(shp.endpoints.uids == ie);
    
    is_cp = find(shp.endpoints.connected_points(is_ep).ids == seq(2));
    ie_cp = find(shp.endpoints.connected_points(ie_ep).ids == seq(end-1));
    
    is_crossp1 = is_cp;
    if is_cp == 1
        is_crossp2 = length(shp.endpoints.connected_points(is_ep).crossp);
    else
        is_crossp2 = is_cp-1;
    end
    
    ie_crossp1 = ie_cp;
    if ie_cp == 1
        ie_crossp2 = length(shp.endpoints.connected_points(ie_ep).crossp);
    else
        ie_crossp2 = ie_cp-1;
    end
    
    crossps1.x = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).x;
    crossps1.y = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).y;
    crossps1.gid = shp.endpoints.connected_points(is_ep).crossp(is_crossp1).gid;
    crossps2.x = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).x;
    crossps2.y = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).y;
    crossps2.gid = shp.endpoints.connected_points(is_ep).crossp(is_crossp2).gid;
    isendpoints = length(shp.endpoints.connected_points(is_ep).crossp) == 2;
    
    crosspe1.x = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).x;
    crosspe1.y = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).y;
    crosspe1.gid = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).gid;
    crosspe2.x = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).x;
    crosspe2.y = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).y;
    crosspe2.gid = shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).gid;
    isendpointe = length(shp.endpoints.connected_points(ie_ep).crossp) == 2;
    
    clnodex = shp.seq(k).cl.coord(:,1)';
    clnodey = shp.seq(k).cl.coord(:,2)';
    clnodewd = shp.seq(k).cl.width;
    
    % Create bank nodes
    nclnodes = length(clnodex);
    newnodes.x_cart = [];
    newnodes.y_cart = [];
    for i=1:nclnodes
        if i == 1
            i1 = 1;
            i2 = 1;
            i3 = 2;
        elseif i == nclnodes
            i1 = nclnodes-1;
            i2 = nclnodes;
            i3 = nclnodes;
        else
            i1 = i-1;
            i2 = i;
            i3 = i+1;
        end
        
        cx = clnodex(i2);
        cy = clnodey(i2);
        
        wd = clnodewd(i2);
        
        if i > 1 && i < length(clnodex)
            x1 = clnodex(i1);
            y1 = clnodey(i1);
            x2 = clnodex(i3);
            y2 = clnodey(i3);
            dx1 = cx-x1;
            dy1 = cy-y1;
            dx2 = x2-cx;
            dy2 = y2-cy;
            costheta = (dx1*dx2+dy1*dy2)/(sqrt(dx1*dx1+dy1*dy1)*sqrt(dx2*dx2+dy2*dy2));
        else
            costheta = 1.0;
        end
        
        coshalftheta = sqrt((1+costheta)/2);
        if coshalftheta < 0.1
            coshalftheta = 0.1;
        end
        
        tvec.x = clnodex(i3) - clnodex(i1);
        tvec.y = clnodey(i3) - clnodey(i1);
        veclen = sqrt(tvec.x^2 + tvec.y^2);
        nvec.x = -tvec.y / veclen * wd * 0.5 / coshalftheta;
        nvec.y =  tvec.x / veclen * wd * 0.5 / coshalftheta;
        newnodes.x_cart(end+1:end+2) = [cx + nvec.x cx - nvec.x];
        newnodes.y_cart(end+1:end+2) = [cy + nvec.y cy - nvec.y];
    end

    % Update crossp coordinates
    if ~isendpoints
        newnodes.x_cart(1:2) = [crossps1.x crossps2.x];
        newnodes.y_cart(1:2) = [crossps1.y crossps2.y];
    else
        shp.endpoints.connected_points(is_ep).crossp(is_crossp1).x = newnodes.x_cart(1);
        shp.endpoints.connected_points(is_ep).crossp(is_crossp2).x = newnodes.x_cart(2);
        shp.endpoints.connected_points(is_ep).crossp(is_crossp1).y = newnodes.y_cart(1);
        shp.endpoints.connected_points(is_ep).crossp(is_crossp2).y = newnodes.y_cart(2);
    end
    if ~isendpointe
        newnodes.x_cart(end-1:end) = [crosspe2.x crosspe1.x];
        newnodes.y_cart(end-1:end) = [crosspe2.y crosspe1.y];
    else
        shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).x = newnodes.x_cart(end);
        shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).x = newnodes.x_cart(end-1);
        shp.endpoints.connected_points(ie_ep).crossp(ie_crossp1).y = newnodes.y_cart(end);
        shp.endpoints.connected_points(ie_ep).crossp(ie_crossp2).y = newnodes.y_cart(end-1);
    end

    % Remove intersecting edges from both ends
    npt = length(newnodes.x_cart);
    if npt > 8
        C = [1 2;(2:2:npt-2)' (4:2:npt)';npt npt-1;flip(npt-1:-2:3)' flip(npt-3:-2:1)'];
        intersects = zeros(1,length(newnodes.x_cart));
        targets = [1,2,size(C,1)/2,size(C,1)/2+1,size(C,1)/2+2,size(C,1)];
        for j=targets
            x1 = newnodes.x_cart(C(j,1));
            y1 = newnodes.y_cart(C(j,1));
            x2 = newnodes.x_cart(C(j,2));
            y2 = newnodes.y_cart(C(j,2));
            for l=1:size(C,1)
                if l==j
                    continue
                end
                x3 = newnodes.x_cart(C(l,1));
                y3 = newnodes.y_cart(C(l,1));
                x4 = newnodes.x_cart(C(l,2));
                y4 = newnodes.y_cart(C(l,2));
                x=[x1 x2 x3 x4];
                y=[y1 y2 y3 y4];
                dt1=det([1,1,1;x(1),x(2),x(3);y(1),y(2),y(3)])*det([1,1,1;x(1),x(2),x(4);y(1),y(2),y(4)]);
                dt2=det([1,1,1;x(1),x(3),x(4);y(1),y(3),y(4)])*det([1,1,1;x(2),x(3),x(4);y(2),y(3),y(4)]);

                if(dt1<-1e-1 && dt2<-1e-1)    %If lines intesect
                    intersects(C(j,:)) = 1;     
                    intersects(C(l,:)) = 1;     
                end
            end
        end
        p_to_del = [];
        if intersects(1) || intersects(2)
            p_to_del = [p_to_del 3 4];
            shp.seq(k).cl.coord(2,:) = [];
            shp.seq(k).cl.depth(2) = [];
            shp.seq(k).cl.width(2) = [];
            shp.seq(k).cl.dists(2) = [];
        end
        if intersects(npt-1) || intersects(npt)
            p_to_del = [p_to_del npt-3 npt-2];
            shp.seq(k).cl.coord(end-1,:) = [];
            shp.seq(k).cl.depth(end-1) = [];
            shp.seq(k).cl.width(end-1) = [];
            shp.seq(k).cl.dists(end-1) = [];
        end
        newnodes.x_cart(p_to_del) = [];
        newnodes.y_cart(p_to_del) = [];
    end

    % Move points that are too close to merger point element
    npt = length(newnodes.x_cart);
    if npt >= 8
        % - end 1
        lxy1 = [newnodes.x_cart(1), newnodes.y_cart(1)];
        lxy2 = [newnodes.x_cart(2), newnodes.y_cart(2)];

        % -- node 1
        pxy = [newnodes.x_cart(3), newnodes.y_cart(3)];
        dist1 = point_to_line_segment_distance(pxy, lxy1, lxy2);

        % -- node 2
        pxy = [newnodes.x_cart(4), newnodes.y_cart(4)];
        dist2 = point_to_line_segment_distance(pxy, lxy1, lxy2);

        if min([dist1, dist2]) < min_spacing
            pxy10 = [newnodes.x_cart(1), newnodes.y_cart(1)];
            pxy12 = [newnodes.x_cart(5), newnodes.y_cart(5)];
            pxy20 = [newnodes.x_cart(2), newnodes.y_cart(2)];
            pxy22 = [newnodes.x_cart(6), newnodes.y_cart(6)];
            newnodes.x_cart(3) = 0.5 * (pxy10(1) + pxy12(1));
            newnodes.y_cart(3) = 0.5 * (pxy10(2) + pxy12(2));
            newnodes.x_cart(4) = 0.5 * (pxy20(1) + pxy22(1));
            newnodes.y_cart(4) = 0.5 * (pxy20(2) + pxy22(2));
        end

        % - end 2
        lxy1 = [newnodes.x_cart(npt-1), newnodes.y_cart(npt-1)];
        lxy2 = [newnodes.x_cart(npt), newnodes.y_cart(npt)];

        % -- node 1
        pxy = [newnodes.x_cart(npt-3), newnodes.y_cart(npt-3)];
        dist1 = point_to_line_segment_distance(pxy, lxy1, lxy2);
        % -- node 2
        pxy = [newnodes.x_cart(npt-2), newnodes.y_cart(npt-2)];
        dist2 = point_to_line_segment_distance(pxy, lxy1, lxy2);

        if min([dist1, dist2]) < min_spacing
            pxy10 = [newnodes.x_cart(npt-1), newnodes.y_cart(npt-1)];
            pxy12 = [newnodes.x_cart(npt-5), newnodes.y_cart(npt-5)];
            pxy20 = [newnodes.x_cart(npt), newnodes.y_cart(npt)];
            pxy22 = [newnodes.x_cart(npt-4), newnodes.y_cart(npt-4)];
            newnodes.x_cart(npt-3) = 0.5 * (pxy10(1) + pxy12(1));
            newnodes.y_cart(npt-3) = 0.5 * (pxy10(2) + pxy12(2));
            newnodes.x_cart(npt-2) = 0.5 * (pxy20(1) + pxy22(1));
            newnodes.y_cart(npt-2) = 0.5 * (pxy20(2) + pxy22(2));
        end
    end

    % Remove points that are too close to merger point element
    npt = length(newnodes.x_cart);
    if npt >= 6
        p_to_del = [];

        % - end 1
        lxy1 = [newnodes.x_cart(1), newnodes.y_cart(1)];
        lxy2 = [newnodes.x_cart(2), newnodes.y_cart(2)];
        llen = norm(lxy1 - lxy2);

        % -- node 1
        pxy = [newnodes.x_cart(3), newnodes.y_cart(3)];
        dist1 = point_to_line_segment_distance(pxy, lxy1, lxy2);

        % -- node 2
        pxy = [newnodes.x_cart(4), newnodes.y_cart(4)];
        dist2 = point_to_line_segment_distance(pxy, lxy1, lxy2);

        % if min([dist1, dist2]) < llen * 0.01
        if min([dist1, dist2]) < min_spacing
            p_to_del = [p_to_del 3 4];
            shp.seq(k).cl.coord(2,:) = [];
            shp.seq(k).cl.depth(2) = [];
            shp.seq(k).cl.width(2) = [];
            shp.seq(k).cl.dists(2) = [];
        end

        % - end 2
        lxy1 = [newnodes.x_cart(npt-1), newnodes.y_cart(npt-1)];
        lxy2 = [newnodes.x_cart(npt), newnodes.y_cart(npt)];
        llen = norm(lxy1 - lxy2);

        % -- node 1
        pxy = [newnodes.x_cart(npt-3), newnodes.y_cart(npt-3)];
        dist1 = point_to_line_segment_distance(pxy, lxy1, lxy2);
        % -- node 2
        pxy = [newnodes.x_cart(npt-2), newnodes.y_cart(npt-2)];
        dist2 = point_to_line_segment_distance(pxy, lxy1, lxy2);

        % if min([dist1, dist2]) < llen * 0.3
        if min([dist1, dist2]) < min_spacing
            p_to_del = [p_to_del npt-3 npt-2];
            shp.seq(k).cl.coord(end-1,:) = [];
            shp.seq(k).cl.depth(end-1) = [];
            shp.seq(k).cl.width(end-1) = [];
            shp.seq(k).cl.dists(end-1) = [];
            removed = true;
        end

        newnodes.x_cart(p_to_del) = [];
        newnodes.y_cart(p_to_del) = [];
    end

    % % Remove tangling edges by distance
    % for ii = 1:3 % Do this correction three times
    %     nlen = length(newnodes.x_cart);
    %     if nlen >= 8
    %         % check distance
    %         p_to_del = [];
    %         % - end 1
    %         p11 = 1; p12 = 3;
    %         x1 = newnodes.x_cart(p11); x2 = newnodes.x_cart(p12);
    %         y1 = newnodes.y_cart(p11); y2 = newnodes.y_cart(p12);
    %         dist1 = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    %         p21 = 2; p22 = 4;
    %         x1 = newnodes.x_cart(p21); x2 = newnodes.x_cart(p22);
    %         y1 = newnodes.y_cart(p21); y2 = newnodes.y_cart(p22);
    %         dist2 = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    %         % if min(dist1,dist2) < shp.seq(k).cl.dists(1)*0.5
    %         if min(dist1,dist2) < min_spacing
    %             del1 = 1;
    %             p_to_del([end+1,end+2]) = [p12,p22];
    %         else
    %             del1 = 0;
    %         end
    % 
    %         % - end 2
    %         p11 = nlen; p12 = nlen-2;
    %         x1 = newnodes.x_cart(p11); x2 = newnodes.x_cart(p12);
    %         y1 = newnodes.y_cart(p11); y2 = newnodes.y_cart(p12);
    %         dist1 = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    %         p21 = nlen-1; p22 = nlen-3;
    %         x1 = newnodes.x_cart(p21); x2 = newnodes.x_cart(p22);
    %         y1 = newnodes.y_cart(p21); y2 = newnodes.y_cart(p22);
    %         dist2 = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    %         % if min(dist1,dist2) < shp.seq(k).cl.dists(end)*0.5
    %         if min(dist1,dist2) < min_spacing
    %             del2 = 1;
    %             p_to_del([end+1,end+2]) = [p12,p22];
    %         else
    %             del2 = 0;
    %         end
    % 
    %         newnodes.x_cart(p_to_del) = [];
    %         newnodes.y_cart(p_to_del) = [];
    % 
    %         if del1
    %             shp.seq(k).cl.coord(2,:) = [];
    %             shp.seq(k).cl.depth(2) = [];
    %             shp.seq(k).cl.width(2) = [];
    %             shp.seq(k).cl.dists(2) = [];
    %         end
    %         if del2
    %             shp.seq(k).cl.coord(end-1,:) = [];
    %             shp.seq(k).cl.depth(end-1) = [];
    %             shp.seq(k).cl.width(end-1) = [];
    %             shp.seq(k).cl.dists(end-1) = [];
    %         end
    %     end
    % end

    [newnodes.y,newnodes.x] = cart2deg(newnodes.x_cart',newnodes.y_cart',lat00,lon00);
    
    depth = [shp.seq(k).cl.depth;shp.seq(k).cl.depth];
    depth = depth(:);
    
    shp.seq(k).bn.coord = [newnodes.x_cart' newnodes.y_cart'];
    shp.seq(k).bn.lonlat = [newnodes.x newnodes.y];
    shp.seq(k).bn.depth = depth;
    shp.seq(k).bn.gid = [crossps1.gid; ...
        crossps2.gid; ...
        (gid_cnt:(gid_cnt+length(newnodes.x_cart)-5))'; ...
        crosspe2.gid; ...
        crosspe1.gid];
    gid_cnt = gid_cnt+length(newnodes.x_cart)-4;
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(S)
        plot(S(i).X_cart,S(i).Y_cart,'k-')
    end
    for i=1:shp.nseq
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
        plot(shp.seq(i).cl.coord(:,1),shp.seq(i).cl.coord(:,2),'*')
        plot(shp.seq(i).bn.coord(:,1),shp.seq(i).bn.coord(:,2),'s')
        for j=1:size(shp.seq(i).cl.coord,1)
            text(shp.seq(i).cl.coord(j,1), shp.seq(i).cl.coord(j,2), int2str(j))
        end
        text(mean(shp.seq(i).bn.coord(:,1)),mean(shp.seq(i).bn.coord(:,2)),int2str(i),color='b')
    end
    axis equal
end

% Add nodes inside large channels
for k=1:shp.nseq
    % For channels
    coord = shp.seq(k).bn.coord;
    clear newnodes
    newnodes.x_cart = [];
    newnodes.y_cart = [];
    newnodes.depth = [];
    shp.seq(k).cl.splitted = zeros(1, size(shp.seq(k).cl.coord,1));
    for i=2:size(coord,1)/2-1
        xy1 = coord(i*2,:);
        xy2 = coord(i*2-1,:);
        dxy = xy2 - xy1;
        wid = norm(dxy);
        nn_add = floor(wid / shp.seq(k).cl.dists(i) - 0.25);
        if nn_add >= 1
            for j=1:nn_add
                xnew = xy1(1) + dxy(1) * j / (nn_add + 1);
                ynew = xy1(2) + dxy(2) * j / (nn_add + 1);
                newnodes.x_cart(end+1) = xnew;
                newnodes.y_cart(end+1) = ynew;
                newnodes.depth(end+1) = shp.seq(k).bn.depth(i*2-1);
            end
            shp.seq(k).cl.splitted(i) = 1;
        end
    end
    if ~isempty(newnodes.x_cart)
        [newnodes.y,newnodes.x] = cart2deg(newnodes.x_cart',newnodes.y_cart',lat00,lon00);
        shp.seq(k).bn.mid_coord = [newnodes.x_cart' newnodes.y_cart'];
        shp.seq(k).bn.mid_lonlat = [newnodes.x newnodes.y];
        shp.seq(k).bn.mid_depth = newnodes.depth';
        shp.seq(k).bn.mid_gid = (gid_cnt:(gid_cnt+length(newnodes.x_cart)-1))';
        gid_cnt = gid_cnt + length(newnodes.x_cart);
    else
        shp.seq(k).bn.mid_coord = [];
        shp.seq(k).bn.mid_lonlat = [];
        shp.seq(k).bn.mid_depth = [];
        shp.seq(k).bn.mid_gid = [];
    end
end

% Add nodes at large merger points
for i=1:length(shp.endpoints.connected_points)
    if isfield(shp.endpoints.connected_points(i).crossp,'mesh')
        shp.endpoints.connected_points(i).crossp = ...
            rmfield(shp.endpoints.connected_points(i).crossp,'mesh');
    end
    crossp = shp.endpoints.connected_points(i).crossp;
    nodes_added = [0 0 0];
    dp_max = crossp(1).depth;
    for j=2:length(crossp)
        dp_max = max(dp_max, crossp(j).depth);
    end
    for j=1:length(crossp)
        if length(crossp) == 2 && j == 1
            shp.endpoints.connected_points(i).crossp(j).mid_coord = [];
            shp.endpoints.connected_points(i).crossp(j).mid_depth = [];
            shp.endpoints.connected_points(i).crossp(j).mid_gid = [];
            continue
        end
        clear newnodes
        j1 = j;
        j2 = mod(j,length(crossp)) + 1;
        xy1 = [crossp(j1).x, crossp(j1).y];
        xy2 = [crossp(j2).x, crossp(j2).y];
        gid1 = crossp(j1).gid;
        gid2 = crossp(j2).gid;
        dxy = xy2 - xy1;
        wid = norm(dxy);
        nn_add = floor(wid / crossp(j).spacing - 0.25);

        newnodes.x_cart = [];
        newnodes.y_cart = [];
        newnodes.depth = [];

        if nn_add >= 1
            for l=1:nn_add
                xnew = xy1(1) + dxy(1) * l / (nn_add + 1);
                ynew = xy1(2) + dxy(2) * l / (nn_add + 1);
                newnodes.x_cart(end+1) = xnew;
                newnodes.y_cart(end+1) = ynew;
                newnodes.depth(end+1) = dp_max;
            end
        end

        if ~isempty(newnodes.x_cart)
            shp.endpoints.connected_points(i).crossp(j).mid_coord = [newnodes.x_cart' newnodes.y_cart'];
            shp.endpoints.connected_points(i).crossp(j).mid_depth = newnodes.depth';
            shp.endpoints.connected_points(i).crossp(j).mid_gid = (gid_cnt:(gid_cnt+length(newnodes.x_cart)-1))';
            gid_cnt = gid_cnt + length(newnodes.x_cart);
            if length(crossp) == 2
                is_end_points(gid1) = 0;
                is_end_points(gid2) = 0;
            end
            nodes_added(j) = length(newnodes.x_cart);
        else
            shp.endpoints.connected_points(i).crossp(j).mid_coord = [];
            shp.endpoints.connected_points(i).crossp(j).mid_depth = [];
            shp.endpoints.connected_points(i).crossp(j).mid_gid = [];
        end
    end
    
    shp.endpoints.connected_points(i).reduced_mid_nodes = false;

    % Collapse added nodes and prepare mesh nodes with gids, coords, and depths here
    if sum(nodes_added ~= 0) == 2
        [min_added, imin_added] = min(nodes_added);
        j1 = mod(imin_added-1+1,3)+1;
        j2 = mod(imin_added-1+2,3)+1;
        if nodes_added(imin_added) ~= nodes_added(j1) && ...
           nodes_added(imin_added) ~= nodes_added(j2)
            mid_coord1 = flip(shp.endpoints.connected_points(i).crossp(j1).mid_coord,1);
            mid_coord2 = shp.endpoints.connected_points(i).crossp(j2).mid_coord;

            mid_n1 = size(mid_coord1,1);
            mid_n2 = size(mid_coord2,1);
            target = 1:min(mid_n1, mid_n2);
            n_not_target1 = mid_n1 - length(target);
            n_not_target2 = mid_n2 - length(target);
            
            mid_coord_new = (mid_coord1(target,:) + mid_coord2(target,:))*0.5;
            mid_coord1(target,:) = mid_coord_new;
            mid_coord2(target,:) = mid_coord_new;
            shp.endpoints.connected_points(i).crossp(j1).mid_coord = flip(mid_coord1,1);
            shp.endpoints.connected_points(i).crossp(j2).mid_coord = mid_coord2;
            
            mid_gid1 = flip(shp.endpoints.connected_points(i).crossp(j1).mid_gid,1);
            mid_gid2 = shp.endpoints.connected_points(i).crossp(j2).mid_gid;
            mid_gid3 = shp.endpoints.connected_points(i).crossp(imin_added).mid_gid;
            if mid_gid1(1) < mid_gid2(1)
                if ~isempty(mid_gid3) && mid_gid2(1) < mid_gid3(1)
                    mid_gid3 = mid_gid3 - length(target);
                end
                shp.endpoints.connected_points(i).crossp(imin_added).mid_gid = mid_gid3;

                mid_gid2(target) = mid_gid1(target);
                mid_gid2(end-n_not_target2+1:end) = mid_gid2(end-n_not_target2+1:end) - length(target);
                shp.endpoints.connected_points(i).crossp(j2).mid_gid = mid_gid2;
            else
                if ~isempty(mid_gid3) && mid_gid1(1) < mid_gid3(1)
                    mid_gid3 = mid_gid3 - length(target);
                end
                shp.endpoints.connected_points(i).crossp(imin_added).mid_gid = mid_gid3;

                mid_gid1(target) = mid_gid2(target);
                mid_gid1(end-n_not_target1+1:end) = mid_gid1(end-n_not_target1+1:end) - length(target);
                shp.endpoints.connected_points(i).crossp(j1).mid_gid = flip(mid_gid1);
            end
            gid_cnt = gid_cnt - length(target);
            shp.endpoints.connected_points(i).reduced_mid_nodes = true;
            shp.endpoints.connected_points(i).mesh.coord = ...
                [shp.endpoints.connected_points(i).crossp(imin_added).x ...
                 shp.endpoints.connected_points(i).crossp(imin_added).y; ...
                 shp.endpoints.connected_points(i).crossp(imin_added).mid_coord; ...
                 shp.endpoints.connected_points(i).crossp(j1).x ...
                 shp.endpoints.connected_points(i).crossp(j1).y; ...
                 shp.endpoints.connected_points(i).crossp(j1).mid_coord(1:n_not_target1+1,:); ...
                 shp.endpoints.connected_points(i).crossp(j2).mid_coord(mid_n2-n_not_target2+1:mid_n2,:)];
            shp.endpoints.connected_points(i).mesh.depth = ones(1, size(shp.endpoints.connected_points(i).mesh.coord,1));
            shp.endpoints.connected_points(i).mesh.depth(:) = dp_max;
            shp.endpoints.connected_points(i).mesh.gid = ...
                [shp.endpoints.connected_points(i).crossp(imin_added).gid; ...
                 shp.endpoints.connected_points(i).crossp(imin_added).mid_gid; ...
                 shp.endpoints.connected_points(i).crossp(j1).gid; ...
                 shp.endpoints.connected_points(i).crossp(j1).mid_gid(1:n_not_target1+1); ...
                 shp.endpoints.connected_points(i).crossp(j2).mid_gid(mid_n2-n_not_target2+1:mid_n2)];
        end
    end
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(S)
        plot(S(i).X_cart,S(i).Y_cart,'k-')
    end
    for i=1:shp.nseq
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'+-')
        plot(shp.seq(i).cl.coord(:,1),shp.seq(i).cl.coord(:,2),'*')
        plot(shp.seq(i).bn.coord(:,1),shp.seq(i).bn.coord(:,2),'s')
    end
    for k=1:shp.nseq
        clnodex = shp.seq(k).cl.coord(:,1)';
        clnodey = shp.seq(k).cl.coord(:,2)';
        plot(clnodex,clnodey,'b-')
    end
    for k=1:shp.nseq
        if ~isempty(shp.seq(k).bn.mid_coord)
            nodex = shp.seq(k).bn.mid_coord(:,1);
            nodey = shp.seq(k).bn.mid_coord(:,2);
            plot(nodex,nodey,'bo')
        end
    end
    for i=1:length(shp.endpoints.connected_points)
        crossp = shp.endpoints.connected_points(i).crossp;
        for j=1:length(crossp)
            nodex = crossp(j).x;
            nodey = crossp(j).y;
            plot(nodex,nodey,'g*')
            % text(nodex,nodey,int2str(j))
            
            if ~isempty(crossp(j).mid_coord)
                mnodex = crossp(j).mid_coord(:,1);
                mnodey = crossp(j).mid_coord(:,2);
                plot(mnodex,mnodey,'ro')
                plot([nodex, mnodex(1)], [nodey, mnodey(1)], '-', color='green', LineWidth=2)
                if length(crossp) >= 3
                    % text(mnodex(1), mnodey(1), int2str(shp.endpoints.connected_points(i).seq(j)))
                end
            end
            for k=1:length(crossp(j).mid_gid)
                text(crossp(j).mid_coord(k,1), crossp(j).mid_coord(k,2), int2str(crossp(j).mid_gid(k)),color='r')
            end
        end
        text(mean([crossp(1).x, crossp(2).x]), mean([crossp(1).y, crossp(2).y]), int2str(i), color='blue')
    end
    for i=1:length(shp.endpoints.uids)
        i0 = shp.endpoints.uids(i);
        x0 = shp.coord(i0,1);
        y0 = shp.coord(i0,2);
        
        if length(shp.endpoints.connected_points(i).ids) == 1
            % Do nothing here
        else
            for j=1:length(shp.endpoints.connected_points(i).ids)
                i1 = shp.endpoints.connected_points(i).ids(j);
                vec1 = shp.endpoints.connected_points(i).vecs(j,:);
                x1 = x0 + vec1(1)*100;
                y1 = y0 + vec1(2)*100;
                plot([x0 x1],[y0 y1],'r-')
            end
        end
    end
    axis equal
end

tPrev = toc_disp(tPrev,'9');

% Generate merger point meshes
clear opts
opts.ref1 = 'preserve';
opts.dtri = 'constrained';
opts.disp = inf;

for i=1:length(shp.endpoints.connected_points)
    if isfield(shp.endpoints.connected_points(i).crossp,'mesh')
        shp.endpoints.connected_points(i).crossp = ...
            rmfield(shp.endpoints.connected_points(i).crossp,'mesh');
    end
    crossp = shp.endpoints.connected_points(i).crossp;
    if length(crossp) < 3
        continue;
    end
    xs = [];
    ys = [];
    dps = [];
    gid = [];
    if shp.endpoints.connected_points(i).reduced_mid_nodes
        assert(isfield(shp.endpoints.connected_points(i), 'mesh'))
        xs = shp.endpoints.connected_points(i).mesh.coord(:,1)';
        ys = shp.endpoints.connected_points(i).mesh.coord(:,2)';
        dps = shp.endpoints.connected_points(i).mesh.depth';
        gid = shp.endpoints.connected_points(i).mesh.gid';
        dps(:) = max(dps);
    else
        for j=1:length(crossp)
            if isempty(crossp(j).mid_coord)
                xs = [xs crossp(j).x];
                ys = [ys crossp(j).y];
                dps = [dps crossp(j).depth];
                gid = [gid crossp(j).gid];
            else
                xs = [xs crossp(j).x crossp(j).mid_coord(:,1)'];
                ys = [ys crossp(j).y crossp(j).mid_coord(:,2)'];
                dps = [dps crossp(j).depth crossp(j).mid_depth'];
                gid = [gid crossp(j).gid crossp(j).mid_gid'];
            end
        end
        dps(:) = max(dps);
    end

    if length(xs) > length(crossp)
        nnd = length(xs);
        edges = [(1:nnd)' [(2:nnd) 1]'];
        [xys,~,elem,~] ...
            = refine2_custom([xs' ys'],edges,[],opts);
        if size(xys,1)-nnd > 0
            ndiff = size(xys,1)-nnd;
            dps = [dps ones(1,ndiff)*max(dps)];
            gid = [gid ((1:ndiff)+gid_cnt-1)];
            gid_cnt = gid_cnt + ndiff;
        end
        xs = xys(:,1)'; ys = xys(:,2)';
        shp.endpoints.connected_points(i).mesh.coord = xys;
        shp.endpoints.connected_points(i).mesh.depth = dps;
        shp.endpoints.connected_points(i).mesh.gid = gid;
    else
        DT = delaunayTriangulation(xs',ys');
        elem = DT.ConnectivityList;
        shp.endpoints.connected_points(i).mesh.coord = [xs' ys'];
        shp.endpoints.connected_points(i).mesh.depth = dps;
        shp.endpoints.connected_points(i).mesh.gid = gid;
    end

    % Remove the tiny elements
    triareas = triangle_area(xs(elem),ys(elem));
    shp.endpoints.connected_points(i).mesh.elem = elem(triareas > 1e-2,:);
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:shp.nseq
        plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'-')
    end
    for i=1:length(shp.endpoints.connected_points)
        crossp = shp.endpoints.connected_points(i);
        if isfield(crossp,'mesh')
            mesh = crossp.mesh;
            if ~isempty(mesh)
                triplot(mesh.elem,mesh.coord(:,1),mesh.coord(:,2))
            end
        end
    end
    for i=1:length(shp.endpoints.uids)
        i0 = shp.endpoints.uids(i);
        x0 = shp.coord(i0,1);
        y0 = shp.coord(i0,2);
        
        if length(shp.endpoints.connected_points(i).ids) == 1
            % Do nothing here
        else
            for j=1:length(shp.endpoints.connected_points(i).ids)
                i1 = shp.endpoints.connected_points(i).ids(j);
                vec1 = shp.endpoints.connected_points(i).vecs(j,:);
                x1 = x0 + vec1(1)*100;
                y1 = y0 + vec1(2)*100;
                plot([x0 x1],[y0 y1],'r-')
            end
        end
    end
    axis equal
end

% Plot for an intermediate check
if plot_level >= 2
        figure;
        hold on;
    %     for i=1:shp.nseq
    %         plot(shp.coord(shp.seq(i).seq,1),shp.coord(shp.seq(i).seq,2),'-')
    %     end
    %     for i=1:length(shp.endpoints.connected_points)
    %         crossp = shp.endpoints.connected_points(i);
    %         if isfield(crossp,'mesh')
    %             mesh = crossp.mesh;
    %             if ~isempty(mesh)
    %                 triplot(mesh.elem,mesh.coord(:,1),mesh.coord(:,2))
    %             end
    %         end
    %     end
        for i=1:shp.nseq
            coord = shp.seq(i).bn.coord;
            npt = size(shp.seq(i).bn.coord,1);
            C = [1 2;(2:2:npt-2)' (4:2:npt)';npt npt-1;flip(npt-1:-2:3)' flip(npt-3:-2:1)'];
            for j=1:size(C,1)
                intrsct = false;
                x1 = coord(C(j,1),1);
                y1 = coord(C(j,1),2);
                x2 = coord(C(j,2),1);
                y2 = coord(C(j,2),2);
                for k=1:size(C,1)
                    if k==j
                        continue
                    end
                    x3 = coord(C(k,1),1);
                    y3 = coord(C(k,1),2);
                    x4 = coord(C(k,2),1);
                    y4 = coord(C(k,2),2);
                    x=[x1 x2 x3 x4];
                    y=[y1 y2 y3 y4];
                    dt1=det([1,1,1;x(1),x(2),x(3);y(1),y(2),y(3)])*det([1,1,1;x(1),x(2),x(4);y(1),y(2),y(4)]);
                    dt2=det([1,1,1;x(1),x(3),x(4);y(1),y(3),y(4)])*det([1,1,1;x(2),x(3),x(4);y(2),y(3),y(4)]);
                    
                    if(dt1<-1e-1 && dt2<-1e-1)
                        intrsct=1;         %If lines intesect
                    end
                end
                if intrsct
                    mkr = 'r-';
                    wdt = 10;
                else
                    mkr = 'k-';
                    wdt = 0.5;
                end
                plot([coord(C(j,1),1),coord(C(j,2),1)],[coord(C(j,1),2),coord(C(j,2),2)],mkr,linewidth=wdt)
            end
            % text(mean([x1,x2]),mean([y1,y2]),int2str(i))
        end
        axis equal
end

tPrev = toc_disp(tPrev,'10');

% Generate channel meshes
clear poly;
for i=1:shp.nseq
    is_cp = shp.seq(i).iconnected_points(1);
    ie_cp = shp.seq(i).iconnected_points(2);
    is_crossp = shp.seq(i).icrossps(1);
    ie_crossp = shp.seq(i).icrossps(2);
    npt_s = size(shp.endpoints.connected_points(is_cp).crossp(is_crossp).mid_coord,1);
    npt_e = size(shp.endpoints.connected_points(ie_cp).crossp(ie_crossp).mid_coord,1);
    npt = size(shp.seq(i).bn.coord,1);
    C = [(1:(npt_s+2)-1)' (2:(npt_s+2))'; ...
         ((npt_s+2):2:npt_s+npt-4)' (npt_s+4:2:npt_s+npt-2)'; ...
         npt_s+npt-2 npt_s+npt+npt_e; ...
         (npt_s+npt+npt_e:-1:npt_s+npt)' (npt_s+npt+npt_e-1:-1:npt_s+npt-1)'; ...
         (npt_s+npt-1:-2:npt_s+5)' (npt_s+npt-3:-2:npt_s+3)';
         npt_s+3 1];
    coord = [shp.seq(i).bn.coord(1,:); ...
             flip(shp.endpoints.connected_points(is_cp).crossp(is_crossp).mid_coord,1); ...
             shp.seq(i).bn.coord(2,:); ...
             shp.seq(i).bn.coord(3:end-1,:); ...
             shp.endpoints.connected_points(ie_cp).crossp(ie_crossp).mid_coord; ...
             shp.seq(i).bn.coord(end,:); ...
             shp.seq(i).bn.mid_coord];
    xs = coord(:,1);
    ys = coord(:,2);

    % if i == 5
    %     figure
    %     hold on
    %     for j=1:size(C,1)
    %         plot(coord(C(j,:),1),coord(C(j,:),2),'-')
    %     end
    %     for j=1:size(coord,1)
    %         text(coord(j,1),coord(j,2),int2str(j))
    %     end
    % end

    DT = delaunayTriangulation(coord(:,1),coord(:,2),C);
    assert(size(DT.Points,1) == size(coord,1))

    P = [C(:,1)', 1];
    polyin = polyshape(coord(P,1)', coord(P,2)', Simplify=false, KeepCollinearPoints=true);
    CC = incenter(DT);
    TF = isinterior(polyin, CC);
    % figure;triplot(DT.ConnectivityList(TF,:), coord(:,1), coord(:,2));axis equal
    
    depth = [shp.seq(i).bn.depth(1); ...
             flip(shp.endpoints.connected_points(is_cp).crossp(is_crossp).mid_depth); ...
             shp.seq(i).bn.depth(2); ...
             shp.seq(i).bn.depth(3:end-1); ...
             shp.endpoints.connected_points(ie_cp).crossp(ie_crossp).mid_depth; ...
             shp.seq(i).bn.depth(end,:); ...
             shp.seq(i).bn.mid_depth];

    gid   = [shp.seq(i).bn.gid(1); ...
             flip(shp.endpoints.connected_points(is_cp).crossp(is_crossp).mid_gid); ...
             shp.seq(i).bn.gid(2); ...
             shp.seq(i).bn.gid(3:end-1); ...
             shp.endpoints.connected_points(ie_cp).crossp(ie_crossp).mid_gid; ...
             shp.seq(i).bn.gid(end,:); ...
             shp.seq(i).bn.mid_gid];

    shp.seq(i).mesh.coord = coord;
    shp.seq(i).mesh.depth = depth';
    shp.seq(i).mesh.gid   = gid';
    shp.seq(i).mesh.boundary = C;
    shp.seq(i).mesh.elem = DT.ConnectivityList(TF, :);

    % Remove tiny elements
    triareas = triangle_area(xs(shp.seq(i).mesh.elem),ys(shp.seq(i).mesh.elem));
    shp.seq(i).mesh.elem = shp.seq(i).mesh.elem(triareas > 1e-2,:);
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:shp.nseq
        mesh = shp.seq(i).mesh;
        triplot(mesh.elem,mesh.coord(:,1),mesh.coord(:,2))
        x = mesh.coord(:,1);
        y = mesh.coord(:,2);
        plot(mean(x(mesh.elem),2),mean(y(mesh.elem),2),'.',Color='r')
        text(mean(x),mean(y),int2str(i))
    end
    for i=1:length(shp.endpoints.connected_points)
        cp = shp.endpoints.connected_points(i);
        if isfield(cp,'mesh')
            mesh = cp.mesh;
            if ~isempty(mesh)
                triplot(mesh.elem,mesh.coord(:,1),mesh.coord(:,2),color='g')
            end
        end
    end
    axis equal
end

% Merge the merger point meshes and the channel meshes

% - Find the maximum gid
max_gid = 0;
for i=1:length(shp.endpoints.connected_points)
    cp = shp.endpoints.connected_points(i);
    if isfield(cp,'mesh') && ~isempty(cp.mesh)
        max_gid = max([max_gid cp.mesh.gid]);
    else
        for j=1:length(cp.crossp)
            max_gid = max([max_gid cp.crossp(j).mid_gid']);
        end
    end
end
for i=1:shp.nseq
    mesh = shp.seq(i).mesh;
    max_gid = max([max_gid mesh.gid]);
end

lonlat = zeros(max_gid, 2);
xy = zeros(max_gid, 2);
depth = ones(1, max_gid) * -99999.0;
elem = [];

% - Merge the merger point meshes
for i=1:length(shp.endpoints.connected_points)
    cp = shp.endpoints.connected_points(i);
    if isfield(cp,'mesh') && ~isempty(cp.mesh)
        mesh_gid = cp.mesh.gid;
        [lat,lon] = cart2deg(cp.mesh.coord(:,1),cp.mesh.coord(:,2),lat00,lon00);
        lonlat(mesh_gid,:) = [lon lat];
        xy(mesh_gid,:) = cp.mesh.coord;
        elem = [elem; mesh_gid(cp.mesh.elem)];
        depth(mesh_gid) = max(cp.mesh.depth);
    else
        gids = [];
        mid_gids = [];
        for j=1:length(cp)
            crossp = cp(j).crossp;
            for k=1:length(crossp)
                [lat,lon] = cart2deg(crossp(k).x,crossp(k).y,lat00,lon00);
                lonlat(crossp(k).gid,:) = [lon lat];
                xy(crossp(k).gid,:) = [crossp(k).x,crossp(k).y];
                depth(crossp(k).gid) = crossp(k).depth;
                gids = [gids crossp(k).gid];
                depth(crossp(k).mid_gid) = crossp(k).mid_depth;
                mid_gids = [mid_gids crossp(k).mid_gid'];
            end
        end
        depth([gids mid_gids]) = max(depth([gids mid_gids]));
    end
end

% - Merge the channel meshes
for i=1:shp.nseq
    mesh = shp.seq(i).mesh;
    [lat,lon] = cart2deg(mesh.coord(:,1),mesh.coord(:,2),lat00,lon00);
    lonlat(mesh.gid,:) = [lon lat];
    xy(mesh.gid,:) = mesh.coord;
    notset = depth(mesh.gid) == -99999.0;
    depth(mesh.gid(notset)) = mesh.depth(notset);
    elem = [elem; mesh.gid(mesh.elem)];
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    patch('faces',elem,'vertices',xy , ...
        'facevertexcdata' , depth', ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    axis equal
end

% Plot for an intermediate check %%%%
if plot_level >= 2
    figure;
    hold on;
    triplot(elem,xy(:,1),xy(:,2));
end

% Plot for an intermediate check %%%%
if plot_level >= 2
    figure;
    hold on;
    triplot(elem,xy(:,1),xy(:,2));
    for i = 1:shp.nseq
        plot(shp.seq(i).cl.coord(:,1),shp.seq(i).cl.coord(:,2),'k-')
    end
    for i=1:length(shp.endpoints.connected_points)
        cp = shp.endpoints.connected_points(i);
        if isfield(cp,'mesh')
            mesh = cp.mesh;
            if ~isempty(mesh)
                triplot(mesh.elem,mesh.coord(:,1),mesh.coord(:,2),'c-')
            end
            for j=1:size(mesh.coord,1)
                text(mesh.coord(j,1),mesh.coord(j,2),int2str(mesh.gid(j)))
            end
        end
    end
    for i=1:length(shp.endpoints.uids)
        i0 = shp.endpoints.uids(i);
        x0 = shp.coord(i0,1);
        y0 = shp.coord(i0,2);
        
        if length(shp.endpoints.connected_points(i).ids) == 1
            % Do nothing here
        else
            for j=1:length(shp.endpoints.connected_points(i).ids)
                i1 = shp.endpoints.connected_points(i).ids(j);
                vec1 = shp.endpoints.connected_points(i).vecs(j,:);
                x1 = x0 + vec1(1)*100;
                y1 = y0 + vec1(2)*100;
                plot([x0 x1],[y0 y1],'r-')
            end
        end
    end
    axis equal
end

tPrev = toc_disp(tPrev,'11-1');

% Make n2e
nnch = size(lonlat,1);
nech = size(elem,1);
n2ech = cell(nnch,1);
for j = 1:nech
    for nm = 1:3
        n2ech{elem(j,nm),1} = [n2ech{elem(j,nm),1} j];
    end
end

% Save the lon-lat at the end points in the channel mesh for later use
endpoints1 = find(is_end_points > 0);
endpoints2 = is_end_points(endpoints1);
mch_lonlat_endpoints1 = lonlat(endpoints1,:);
mch_lonlat_endpoints2 = lonlat(endpoints2,:);

% Collapse endpoint nodes
nodes_toberemoved = [];
elems_toberemoved = [];
lonlat_new = lonlat;
lonlat_collapsed = [];
exyrange = [];
for i = find(is_end_points > 0)
    i1 = i;
    i2 = is_end_points(i);

    lonp = mean(lonlat([i1,i2],1));
    latp = mean(lonlat([i1,i2],2));
    [ielem,exyrange] = find_element_containing_xy(msub.p,msub.t,lonp,latp,exyrange);
    if ielem <= 0 % do not remove this end if it is outside the orignal grid domain
        continue;
    end

    nelem1 = n2ech{i1};
    nelem2 = n2ech{i2};
    if length(nelem1) ~= 1 && length(nelem2) ~= 1 % do not remove this end because there are more than 2 elements associated.
        continue
    end
    if length(nelem1) ~= 1
        ii1 = i1;
        ii2 = i2;
        elem_toberemoved = nelem2;
    else
        ii1 = i2;
        ii2 = i1;
        elem_toberemoved = nelem1;
    end

    clon = mean(lonlat([i1 i2],1));
    clat = mean(lonlat([i1 i2],2));
    lonlat_new(ii1,1) = clon;
    lonlat_new(ii1,2) = clat;
    lonlat_collapsed(end+1,1) = clon;
    lonlat_collapsed(end,2) = clat;
    nodes_toberemoved(end+1) = ii2;
    elems_toberemoved(end+1) = elem_toberemoved;
end

tPrev = toc_disp(tPrev,'11-2');

% Register the nodes to be removed in the shp struct
for i=1:length(shp.endpoints.connected_points)
    cp = shp.endpoints.connected_points(i);
    for j=1:length(cp.crossp)
        shp.endpoints.connected_points(i).crossp(j).removed ...
            = ismember(shp.endpoints.connected_points(i).crossp(j).gid,nodes_toberemoved);
    end
end

% Remove nodes without lonlat assigned
nodes_toberemoved = [nodes_toberemoved find(lonlat_new(:,1) == 0)];

% Remove nodes and elements
nodes_tostay = ~ismember(1:nnch,nodes_toberemoved);
elems_tostay = ~ismember(1:nech,elems_toberemoved);

lonlat_new = lonlat_new(nodes_tostay,:);
map_node = zeros(size(lonlat,1),1);
map_node(nodes_tostay) = 1:nnz(nodes_tostay);
elem_new = map_node(elem(:,:));
elem_new = elem_new(elems_tostay,:);
depth_new = depth(nodes_tostay);

tPrev = toc_disp(tPrev,'11');

% Update gid
% - Merge point meshes
for i=1:length(shp.endpoints.connected_points)
    cp = shp.endpoints.connected_points(i);
    for j=1:length(cp)
        crossp = cp(j).crossp;
        for k=1:length(crossp)
            shp.endpoints.connected_points(i).crossp(k).gid ...
                = map_node(shp.endpoints.connected_points(i).crossp(k).gid);
            shp.endpoints.connected_points(i).crossp(k).mid_gid ...
                = map_node(shp.endpoints.connected_points(i).crossp(k).mid_gid);
        end
    end
    if isfield(cp,'mesh') && ~isempty(cp.mesh)
        shp.endpoints.connected_points(i).mesh.gid ...
            = map_node(shp.endpoints.connected_points(i).mesh.gid);
    end
end

% - Channel meshes
for i=1:shp.nseq
    shp.seq(i).mesh.gid ...
        = map_node(shp.seq(i).mesh.gid);
end

% Set the msh object
mch = msh();
mch.p = lonlat_new;
mch.b = depth_new';
mch.t = elem_new;
nnch = size(mch.p,1);

% Plot for an intermediate check
if plot_level >= 2
    figure;
    triplot(mch.t,mch.p(:,1),mch.p(:,2));
    hold on
    nid = setdiff((1:size(mch.p,1))', unique(mch.t));
    disp(nid)
    plot(mch.p(nid,1), mch.p(nid,2), 'ro')
    axis equal
end

% Write the channel mesh
mch.write([out14 '_mch'],{'14'})

% Plot for an intermediate check
if plot_level >= 2
    figure;
    patch('faces',mch.t,'vertices',mch.p , ...
        'facevertexcdata' , mch.b, ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    axis equal
end

tPrev = toc_disp(tPrev,'12');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Embed the channel mesh to a larger mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nd_in_stay = (1:nnsub)';
nd_in_removed = [];

lon1_prev = 0.0;
lat1_prev = 0.0;

for j = 1:size(shp.coord,1)
    lon1 = shp.lonlat(j,1);
    lat1 = shp.lonlat(j,2);
    dl = shp.dists(j)*removal_range;
    [x1,y1] = deg2cart(lat1,lon1,lat00,lon00);
    [dum,lond] = cart2deg(x1+dl,y1,lat00,lon00);
    [latd,dum] = cart2deg(x1,y1+dl,lat00,lon00);
    dlon = lond - lon1;
    dlat = latd - lat1;
    plon = msub.p(nd_in_stay,1);
    plat = msub.p(nd_in_stay,2);
    targets = plon >= lon1-dlon & plon <= lon1+dlon & plat >= lat1-dlat & plat <= lat1+dlat;

    d_prev = m_lldist([lon1_prev lon1],[lat1_prev,lat1])*1000.0;
    if d_prev < 0.25*dl
        continue;
    end

    for k = nd_in_stay(targets)'
        lon2 = msub.p(k,1);
        lat2 = msub.p(k,2);
        d = m_lldist([lon1 lon2],[lat1 lat2])*1000.0;
        if d < dl
            nd_in_removed = [nd_in_removed k];
        end
    end
    nd_in_stay = setdiff(nd_in_stay,nd_in_removed);
end

nd_in_removed = nd_in(unique(nd_in_removed));

tPrev = toc_disp(tPrev,'12.3');

% nd_in_removed = setdiff(nd_in_removed,nodes_tobekept); % Eliminate boundary nodes from the removed node list
% el_in_stay_mem1 = ~ismember(m.t(:,1),nd_in_removed);
% el_in_stay_mem2 = ~ismember(m.t(:,2),nd_in_removed);
% el_in_stay_mem3 = ~ismember(m.t(:,3),nd_in_removed);
% el_in_stay_mem = el_in_mem & el_in_stay_mem1 & el_in_stay_mem2 & el_in_stay_mem3;
% el_in_stay = find(el_in_stay_mem);
% nd_in_stay = m.t(el_in_stay,:); nd_in_stay = sort(unique(nd_in_stay(:)));
% nd_in_removed = setdiff(nd_in,nd_in_stay);
%
nd_in_removed = setdiff(nd_in_removed,nodes_tobekept); % Eliminate boundary nodes from the removed node list
nd_in_stay = setdiff(nd_in,nd_in_removed);
el_in_stay_mem1 = ~ismember(m.t(:,1),nd_in_removed);
el_in_stay_mem2 = ~ismember(m.t(:,2),nd_in_removed);
el_in_stay_mem3 = ~ismember(m.t(:,3),nd_in_removed);
el_in_stay_mem = el_in_mem & el_in_stay_mem1 & el_in_stay_mem2 & el_in_stay_mem3;
el_in_stay = find(el_in_stay_mem);

tPrev = toc_disp(tPrev,'12.5');

% Create a subset of the subset mesh
msub_stay = msh();
msub_stay.t = m.t(el_in_stay,:);
nodes_sub_stay = sort(unique(msub_stay.t(:)));

% Plot for an intermediate check
if plot_level >= 2
    figure;
    triplot(msub.t,msub.p(:,1),msub.p(:,2));
    triplot(m.t(el_in_stay,:),m.p(:,1),m.p(:,2));
    hold on
    plot(m.p(nd_in_removed,1),m.p(nd_in_removed,2),'o')
    plot(m.p(nodes_tobekept,1),m.p(nodes_tobekept,2),'+')
    plot(m.p(nodes_sub_stay,1),m.p(nodes_sub_stay,2),'x',LineWidth=2)
    axis equal
end

tPrev = toc_disp(tPrev,'12.6');

for j = 1:size(msub_stay.t,1)
    for nm = 1:3
        newi = find(nodes_sub_stay == msub_stay.t(j,nm));
        msub_stay.t(j,nm) = newi;
    end
end

tPrev = toc_disp(tPrev,'12.7');

msub_stay.p = m.p(nodes_sub_stay,:);
msub_stay.b = m.b(nodes_sub_stay);
nnsub_stay = size(msub_stay.p,1);
nesub_stay = size(msub_stay.t,1);

% Plot for an intermediate check
if plot_level >= 2
    figure
    triplot(msub.t,msub.p(:,1),msub.p(:,2))
    hold on
    triplot(msub_stay.t,msub_stay.p(:,1),msub_stay.p(:,2),'r')
    plot(shp.lonlat(:,1),shp.lonlat(:,2),'o')
    axis equal
end

tPrev = toc_disp(tPrev,'13');

% Find elements along the interior boundary on the removed side
el_in_removed_mem = ~el_in_stay_mem & el_in_mem;
el_in_removed_on_boundary_mem = el_in_removed_mem & ...
    (ismember(m.t(:,1),nd_in_stay) | ...
    ismember(m.t(:,2),nd_in_stay) | ...
    ismember(m.t(:,3),nd_in_stay));
% Find boundary nodes
nd_boundary_el = m.t(el_in_removed_on_boundary_mem,:);
nd_boundary_el = unique(nd_boundary_el(:));
nd_boundary = setdiff(nd_boundary_el,nd_in_removed);
nd_boundary0 = nd_boundary;

% Plot for an intermediate check
if plot_level >= 2
    figure
    triplot(m.t(el_in_removed_on_boundary_mem,:),m.p(:,1),m.p(:,2))
    hold on
    plot(shp.lonlat(:,1),shp.lonlat(:,2),'o')
    plot(m.p(nd_boundary,1),m.p(nd_boundary,2),'go')
    plot(m.p(nd_in_removed,1),m.p(nd_in_removed,2),'o')
    axis equal
end

% Create a sequence of boundary nodes (= polygon)
tsub = m.t(el_in_removed_on_boundary_mem,:);
nd_boundary_poly = nd_boundary(1);
nd_boundary_remaining = [nd_boundary(2:end);nd_boundary(1)];
j = 1;
while true
    if j > length(nd_boundary_poly)
        break;
    end
    i0 = nd_boundary_poly(j);
    found = false;
    for i1 = nd_boundary_remaining'
        el_match = ...
            (tsub(:,1) == i0 & tsub(:,2) == i1) | ...
            (tsub(:,2) == i0 & tsub(:,3) == i1) | ...
            (tsub(:,3) == i0 & tsub(:,1) == i1);
        assert(sum(el_match) <= 1)
        if sum(el_match) == 1
            nd_boundary_poly = [nd_boundary_poly i1];
            found = true;
            break;
        end
    end
    if ~found
        nd_boundary(nd_boundary==nd_boundary_poly(end)) = [];
        nd_boundary_poly(end) = [];
        j = j - 1;
    else
        nd_boundary_remaining(nd_boundary_remaining==i1) = [];
        if nd_boundary_poly(end) == nd_boundary_poly(1)
            break;
        end
        j = j + 1;
    end
end
nd_boundary = setdiff(nd_boundary, nd_boundary_remaining);

% Plot for an intermediate check
if plot_level >= 2
    figure
    triplot(m.t(el_in_removed_on_boundary_mem,:),m.p(:,1),m.p(:,2))
    hold on
    plot(shp.lonlat(:,1),shp.lonlat(:,2),'k.')
    plot(m.p(nd_boundary,1),m.p(nd_boundary,2),'o',LineWidth=2)
    plot(m.p(nd_boundary_poly,1),m.p(nd_boundary_poly,2),'-',LineWidth=2)
    axis equal
end

% Create polylines
nn = size(m.p,1);
m.coord = [];
[m.coord(1:nn,1),m.coord(1:nn,2)] ...
    = deg2cart(m.p(:,2),m.p(:,1),lat00,lon00);
exyrange = [];
polys = {};
polys{1} = m.coord(nd_boundary_poly,:);
while true
    found = false;
    for j=1:length(polys)
        poly = polys{j};
        for i=1:size(poly,1)-1
            x1 = poly(i,1);
            y1 = poly(i,2);
            x2 = poly(i+1,1);
            y2 = poly(i+1,2);
            a1 = (y2-y1) / (x2 - x1);
            b1 = y1 - a1*x1;
            for jseq=1:shp.nseq
                for icl=1:size(shp.seq(jseq).cl.coord,1)-1
                    xcl1 = shp.seq(jseq).cl.coord(icl,1);
                    ycl1 = shp.seq(jseq).cl.coord(icl,2);
                    xcl2 = shp.seq(jseq).cl.coord(icl+1,1);
                    ycl2 = shp.seq(jseq).cl.coord(icl+1,2);
                    a2 = (ycl2-ycl1) / (xcl2 - xcl1);
                    b2 = ycl1 - a2*xcl1;
                    X = (b1-b2)/(a2-a1);
                    cond1 = x1 <= X && X <= x2;
                    cond2 = x1 >= X && X >= x2;
                    cond3 = xcl1 <= X && X <= xcl2;
                    cond4 = xcl1 >= X && X >= xcl2;
                    if ( cond1 || cond2 ) && ( cond3 || cond4 )
                        [elem,exyrange] = find_element_containing_xy(msub.coord,msub.t,xcl1,ycl1,exyrange);
                        if elem > 0
                            ibn1 = (icl-1)*2+1;
                            ibn2 = (icl-1)*2+2;
                            icla = icl;
                            % shp.seq(jseq).indomain(icl+1:end) = 0;
                        else
                            ibn1 = (icl)*2+1;
                            ibn2 = (icl)*2+2;
                            icla = icl+1;
                            % shp.seq(jseq).indomain(1:icl) = 0;
                        end
                        xbn1 = shp.seq(jseq).bn.coord(ibn1,1);
                        ybn1 = shp.seq(jseq).bn.coord(ibn1,2);
                        xbn2 = shp.seq(jseq).bn.coord(ibn2,1);
                        ybn2 = shp.seq(jseq).bn.coord(ibn2,2);
                        nx0 = -(ycl2-ycl1);
                        ny0 = xcl2-xcl1;
                        dx1 = x1-xcl1;
                        dy1 = y1-ycl1;
                        dxbn1 = xbn1-xcl1;
                        dybn1 = ybn1-ycl1;
                        dot1 = dx1*nx0 + dy1*ny0;
                        dotbn1 = dxbn1*nx0 + dybn1*ny0;
                        if dot1 * dotbn1 >= 0
                            x1new = xbn1;
                            y1new = ybn1;
                            x2new = xbn2;
                            y2new = ybn2;
                        else
                            x1new = xbn2;
                            y1new = ybn2;
                            x2new = xbn1;
                            y2new = ybn1;
                        end
                        poly1 = [poly(1:i,:); x1new y1new];
                        poly2 = [x2new y2new; poly(i+1:end,:)];
                        while true
                            dist = sqrt((x1new-poly1(end-1,1))^2 + (y1new-poly1(end-1,2))^2);
                            if dist < shp.seq(jseq).cl.dists(icla) && size(poly1,1) >= 3
                                ii = m.coord(nd_boundary,1) == poly1(end-1,1) & m.coord(nd_boundary,2) == poly1(end-1,2);
                                nd_boundary(ii) = [];
                                poly1(end-1,:) = [];
                            else
                                break
                            end
                        end
                        while true
                            dist = sqrt((x2new-poly2(2,1))^2 + (y2new-poly2(2,2))^2);
                            if dist < shp.seq(jseq).cl.dists(icla) && size(poly2,1) >= 3
                                ii = m.coord(nd_boundary,1) == poly2(2,1) & m.coord(nd_boundary,2) == poly2(2,2);
                                nd_boundary(ii) = [];
                                poly2(2,:) = [];
                            else
                                break
                            end
                        end
                        polys{j} = poly1;
                        polys{end+1} = poly2;
                        
                        found = true;
                    end
                    if found
                        break
                    end
                end
                if found
                    break
                end
            end
            if found
                break
            end
        end
        if found
            break
        end
    end
    if ~found
        break
    end
end

polys1 = {};
polys2 = {};
j = 1;
for jseq=1:shp.nseq
    shp.seq(jseq).indomain = ones(1,size(shp.seq(jseq).cl.coord,1));
    for icl=1:size(shp.seq(jseq).cl.coord,1)
        xcl = shp.seq(jseq).cl.coord(icl,1);
        ycl = shp.seq(jseq).cl.coord(icl,2);
        
        [elem,exyrange] = find_element_containing_xy(msub.coord,msub.t,xcl,ycl,exyrange);

        if elem <= 0
            shp.seq(jseq).indomain(icl) = 0;
        end
    end
    poly1 = zeros(0,2);
    poly2 = zeros(0,2);

    % Add mid nodes
    if shp.seq(jseq).indomain(1) == 1
        icp = shp.seq(jseq).iconnected_points(1);
        if length(shp.endpoints.connected_points(icp).crossp) == 2
            icrossp = shp.seq(jseq).icrossps(1);
            if ~isempty(shp.endpoints.connected_points(icp).crossp(icrossp).mid_coord)
                mid_coord = shp.endpoints.connected_points(icp).crossp(icrossp).mid_coord;
                ibn1 = 1;
                poly2(end+1,:) = shp.seq(jseq).bn.coord(ibn1,:);
                poly2(end+1:end+size(mid_coord,1),:) = flip(mid_coord,1);
            end
        end
    end
    
    for icl=1:size(shp.seq(jseq).cl.coord,1)
        if shp.seq(jseq).indomain(icl) == 1
            bn = shp.seq(jseq).bn.coord;
            ibn1 = (icl-1)*2+1;
            ibn2 = (icl-1)*2+2;
            if ismember(shp.seq(jseq).bn.gid(ibn1),abs(is_end_points))
                poly1(end+1,:) = mean([bn(ibn1,:);bn(ibn2,:)]);
                poly2(end+1,:) = poly1(end,:);
            else
                poly1(end+1,:) = bn(ibn1,:);
                poly2(end+1,:) = bn(ibn2,:);
            end
        else
            if ~isempty(poly1)
                polys1{j} = poly1;
                polys2{j} = flip(poly2,1);
                poly1 = zeros(0,2);
                poly2 = zeros(0,2);
                j = j + 1;
            end
        end
    end

    % Add mid nodes
    if ~isempty(poly1)
        if shp.seq(jseq).indomain(size(shp.seq(jseq).cl.coord,1)) == 1
            icp = shp.seq(jseq).iconnected_points(2);
            if length(shp.endpoints.connected_points(icp).crossp) == 2
                icrossp = shp.seq(jseq).icrossps(2);
                if ~isempty(shp.endpoints.connected_points(icp).crossp(icrossp).mid_coord)
                    mid_coord = shp.endpoints.connected_points(icp).crossp(icrossp).mid_coord;
                    poly1(end+1:end+size(mid_coord,1),:) = mid_coord;
                    poly1(end+1,:) = bn(end,:);
                end
            end
        end
    end

    if size(poly1,1) >= 2
        polys1{j} = poly1;
        polys2{j} = flip(poly2,1);
        j = j + 1;
    end
end

% Plot for an intermediate check
[mch.coord(1:nnch,1),mch.coord(1:nnch,2)] ...
    = deg2cart(mch.p(:,2),mch.p(:,1),lat00,lon00);
if plot_level >= 2
    figure;
    triplot(msub.t,msub.coord(:,1),msub.coord(:,2))
    hold on
    triplot(mch.t,mch.coord(:,1),mch.coord(:,2))
%     for j=1:length(polys)
%         poly = polys{j};
%         plot(poly(:,1),poly(:,2),'o-');
%     end
    for j=1:length(polys1)
        poly = polys1{j};
        plot(poly(:,1),poly(:,2),'ro-',LineWidth=2);
        for k=1:size(poly,1)
            text(poly(k,1),poly(k,2),int2str(k))
        end

        poly = polys2{j};
        plot(poly(:,1),poly(:,2),'co-',LineWidth=2);
        for k=1:size(poly,1)
            text(poly(k,1),poly(k,2),int2str(k))
        end
    end
    axis equal
end

% Outer polylines
poly_outer = m.p(nd_boundary_poly,:);
nn_outer = size(poly_outer,1);

% Set the cartesian coordinates
[msub.coord(1:nnsub,1),msub.coord(1:nnsub,2)] ...
    = deg2cart(msub.p(:,2),msub.p(:,1),lat00,lon00);
[poly_outer_cart(1:nn_outer,1),poly_outer_cart(1:nn_outer,2)] ...
    = deg2cart(poly_outer(:,2),poly_outer(:,1),lat00,lon00);

% Make the list of the representative mesh size at the original mesh nodes
sub_nd_mesh_sizes = zeros(1,nnsub);
for i=1:nnsub
    %     subbarlens1 = subbarlens(subbars(:,1) == i);
    %     subbarlens2 = subbarlens(subbars(:,2) == i);
    %     distmin = min([subbarlens1;subbarlens2]);
    x0 = msub.coord(i,1);
    y0 = msub.coord(i,2);
    
    % Look for the closest channel node
    l_chnd = sqrt((shp.coord(:,1)-x0).^2 + (shp.coord(:,2)-y0).^2);
    [min_l_chnd,min_i_l_chnd] = min(l_chnd);
    min_dist_chnd = shp.dists(min_i_l_chnd);
    
    % Look for the closest outer polyline node
    l_opnd = sqrt((poly_outer_cart(:,1)-x0).^2 + (poly_outer_cart(:,2)-y0).^2);
    [min_l_opnd,min_i_l_opnd] = min(l_opnd);
    min_i_l_opnd_gid = find(msub.p(:,1)==poly_outer(min_i_l_opnd,1) & ...
        msub.p(:,2)==poly_outer(min_i_l_opnd,2));
    if length(min_i_l_opnd_gid) > 1
        min_i_l_opnd_gid = min_i_l_opnd_gid(1);
    end
    subbarlens1 = subbarlens(subbars(:,1) == min_i_l_opnd_gid);
    subbarlens2 = subbarlens(subbars(:,2) == min_i_l_opnd_gid);
    min_dist_opnd = min([subbarlens1;subbarlens2]);
    
    sum_min_l = min_l_chnd + min_l_opnd;
    w = min_l_chnd / sum_min_l;
    distmin = (1 - w)*min_dist_chnd + w*min_dist_opnd;
    
    sub_nd_mesh_sizes(i) = distmin*rep_dist_adjust;
end
% figure
% for i=1:10:nnsub
%     circle(msub.coord(i,1),msub.coord(i,2),sub_nd_mesh_sizes(i));
% end
% axis equal

tPrev = toc_disp(tPrev,'15');

% merge polylines
% - consolidate them into one cell arrays
polys_all = {};
for i=1:length(polys)
    polys_all{end+1} = polys{i};
end
if ~no_channelmesh
    for i=1:length(polys1)
        polys_all{end+1} = polys1{i};
    end
    for i=1:length(polys2)
        polys_all{end+1} = polys2{i};
    end
end

if plot_level >= 2
    figure;
    triplot(msub.t,msub.coord(:,1),msub.coord(:,2))
    hold on
    triplot(mch.t,mch.coord(:,1),mch.coord(:,2))
    for j=1:length(polys_all)
        poly = polys_all{j};
        plot(poly(:,1),poly(:,2),'o-',LineWidth=2);
        % plot(poly(2,1),poly(2,2),'k*',LineWidth=2);
    end
    axis equal
end

% - connecting polyines at heads and tails
while true
    updated = false;
    for i=1:length(polys_all)
        polyi = polys_all{i};
        for j=i+1:length(polys_all)
            polyj = polys_all{j};
            if  polyi(end,1) == polyj(1,1) && polyi(end,2) == polyj(1,2)
                polyi_new = [polyi(1:end-1,:);polyj];
                polys_all{i} = polyi_new;
                polys_all(j) = [];
                updated = true;
                break;
            end
        end
        if updated
            break;
        end
    end
    if ~updated
        break;
    end
end

if plot_level >= 2
    figure;
    triplot(msub.t,msub.coord(:,1),msub.coord(:,2))
    hold on
    triplot(mch.t,mch.coord(:,1),mch.coord(:,2))
    for j=1:length(polys_all)
        poly = polys_all{j};
        plot(poly(:,1),poly(:,2),'o-',LineWidth=2);
        % for k=1:size(poly,1)
        %     text(poly(k,1), poly(k,2), int2str(k))
        % end
    end
    axis equal
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    hold on;
    for i=1:length(polys_all)
        poly = polys_all{i};
        coord = poly;
        npt = size(poly,1);
        C = [(1:npt-1)' (2:npt)'];
        for j=1:size(C,1)
            intrsct = false;
            x1 = coord(C(j,1),1);
            y1 = coord(C(j,1),2);
            x2 = coord(C(j,2),1);
            y2 = coord(C(j,2),2);
            for k=1:size(C,1)
                if k==j
                    continue
                end
                x3 = coord(C(k,1),1);
                y3 = coord(C(k,1),2);
                x4 = coord(C(k,2),1);
                y4 = coord(C(k,2),2);
                x=[x1 x2 x3 x4];
                y=[y1 y2 y3 y4];
                dt1=det([1,1,1;x(1),x(2),x(3);y(1),y(2),y(3)])*det([1,1,1;x(1),x(2),x(4);y(1),y(2),y(4)]);
                dt2=det([1,1,1;x(1),x(3),x(4);y(1),y(3),y(4)])*det([1,1,1;x(2),x(3),x(4);y(2),y(3),y(4)]);
                
                if(dt1<-1e-1 && dt2<-1e-1)
                    intrsct=1;         %If lines intesect
                end
            end
            if intrsct
                mkr = 'r-';
                wdt = 10;
            else
                mkr = 'k-';
                wdt = 0.5;
            end
            plot([coord(C(j,1),1),coord(C(j,2),1)],[coord(C(j,1),2),coord(C(j,2),2)],mkr,linewidth=wdt)
        end
        text(mean([x1,x2]),mean([y1,y2]),int2str(i))
    end
    axis equal
end

% Build the connecting mesh
% mconn = build_mesh_between_polylines(poly_outer, poly_inner, rseed, nnsub, msub, sub_nd_mesh_sizes, reject_margin, max_num_rejected, niter_relax, channel_spacing, msub_nodes_inchannel, neinodessub, ignore_depth_in_channel, no_channelmesh);
mconn = build_mesh_inside_polylines(polys_all, lat00, lon00, rseed, nnsub, msub, sub_nd_mesh_sizes, reject_margin, max_num_rejected, niter_relax, channel_spacing, msub_nodes_inchannel, neinodessub, ignore_depth_in_channel, no_channelmesh, plot_level);

% Update the coordinates to match the original along the outer boundary and
% the channel mesh
for ind = nd_boundary'
    lonlati = m.p(ind,:);

    jnd = dsearchn(mconn.p,double(mconn.t),lonlati);
    lonlatj = mconn.p(jnd,:);

    if m_lldist([lonlati(1) lonlatj(1)],[lonlati(2) lonlatj(2)])*1000.0 < 0.1
        mconn.p(jnd,:) = lonlati;
    end
end
for ind = 1:size(mch.p,1)
    lonlati = mch.p(ind,:);

    jnd = dsearchn(mconn.p,double(mconn.t),lonlati);
    lonlatj = mconn.p(jnd,:);

    if m_lldist([lonlati(1) lonlatj(1)],[lonlati(2) lonlatj(2)])*1000.0 < 0.1
        mconn.p(jnd,:) = lonlati;
    end
end

if ~no_channelmesh
    % Remove the collapsed node elements
    nodes_toberemoved_ch = zeros(1,size(mch.p,1));
    
    for icnd = 1:size(lonlat_collapsed,1)
        lonlati = lonlat_collapsed(icnd,:);
    
        jch = dsearchn(mch.p,double(mch.t),lonlati,[]);
        jconn = dsearchn(mconn.p,double(mconn.t),lonlati,[]);
    
        assert(~isempty(jch))
        assert(~isempty(jconn))
    
        % Remove the element that has node jch as one of its nodes
        nodes_toberemoved_ch(jch) = 1;
        elem_ch = find(mch.t(:,1) == jch | mch.t(:,2) == jch | mch.t(:,3) == jch);
        assert(length(elem_ch) == 1)
        jj = find(mch.t(elem_ch,:) == jch);
    
        tbl1 = [2 3 1];
        tbl2 = [3 1 2];
    
        jj1 = tbl1(jj);
        lonlatjj1 = mch.p(mch.t(elem_ch,jj1),:);
        jconn1 = dsearchn(mconn.p,double(mconn.t),lonlatjj1,[]);
    
        jj2 = tbl2(jj);
        lonlatjj2 = mch.p(mch.t(elem_ch,jj2),:);
        jconn2 = dsearchn(mconn.p,double(mconn.t),lonlatjj2,[]);
    
        new_elem_conn = [jconn, jconn1, jconn2];
        mconn.t(end+1,1:3) = new_elem_conn; % Add a new element to mconn
        mch.t(elem_ch,:) = []; % Remove an emlement from mch
    
        if indomain_flowbc
            mch_lonlat_endpoints1(end+1,:) = lonlatjj1;
            mch_lonlat_endpoints2(end+1,:) = lonlatjj2;
        end
    end
    %- Remove nodes from mch
    map_remove_collapsed = zeros(size(mch.p,1),2);
    map_remove_collapsed(:,1) = 1:size(mch.p,1);
    map_remove_collapsed(nodes_toberemoved_ch==0,2) = 1:(size(mch.p,1)-sum(nodes_toberemoved_ch));
    mch.p(nodes_toberemoved_ch==1,:) = [];
    mch.b(nodes_toberemoved_ch==1,:) = [];
    map_remove_collapsed2 = map_remove_collapsed(:,2);
    mch.t = map_remove_collapsed2(mch.t);
end

if write_output %%%%%%%%%%%%%%%%%%%
    mconn.write('mconn');
end

nnconn = size(mconn.p,1);
neconn = size(mconn.t,1);

tPrev = toc_disp(tPrev,'16');

% Merge the channel grid

% Shift the node IDs of the channel mesh
mch_new = mch;
mch_new.t = mch_new.t + nnconn;

% Put a new set of channel mesh nodes
mch_new.p = mch.p;

% Add the channel mesh nodes to the connecting mesh nodes
if no_channelmesh
    mconnch = mconn;
else
    mconnch = msh();
    mconnch.p = [mconn.p;mch_new.p];
    mconnch.t = [mconn.t;mch_new.t];
    mconnch.b = [mconn.b;mch_new.b];
end

% # of nodes and elements
nnconnch = size(mconnch.p,1);
neconnch = size(mconnch.t,1);

% Plot for an intermediate check
if plot_level >= 2
    figure
    triplot(mconnch.t,mconnch.p(:,1),mconnch.p(:,2))
    axis equal
end

% Create the boundary strings
% - Create the channel mesh boundary polylines
% -- Make n2e
nnch = size(mch.p,1);
nech = size(mch.t,1);
n2ech = cell(nnch,1);
for j = 1:nech
    for nm = 1:3
        n2ech{mch.t(j,nm)} = [n2ech{mch.t(j,nm)} j];
    end
end
for i = 1:nnch
    n2ech{i} = unique(n2ech{i});
end

tPrev = toc_disp(tPrev,'16-1');

% -- Make boundary polylines
i_start = 1;
i = i_start;
j = 1;
polyi = {};
polyi{j} = i; % first node
icheckedout = zeros(size(mch.p,1),1);
icheckedout(i) = 1;
while true
    neiele = n2ech{i};
    found = false;
    for e = neiele
        nm = mch.t(e,:);
        ii = find(nm == i);
        ii_next = mod(ii+1-1,3)+1;
        i_next = nm(ii_next);
        neiele_next = n2ech{i_next};
        if length(intersect(neiele,neiele_next)) == 1
            found = true;
            break;
        end
    end

    assert(found)
    
    polyi{j} = [polyi{j} i_next];
    i = i_next;
    icheckedout(i) = 1;
    if i_next == i_start
        if all(icheckedout)
            break;
        else
            while true
                i = find(icheckedout == 0, 1);
                icheckedout(i) = 1;
                neiele = n2ech{i};
                found = false;
                for e = neiele
                    nm = mch.t(e,:);
                    ii = find(nm == i);
                    ii_next = mod(ii+1-1,3)+1;
                    i_next = nm(ii_next);
                    neiele_next = n2ech{i_next};
                    if length(intersect(neiele,neiele_next)) == 1
                        found = true;
                        break;
                    end
                end
                if found || all(icheckedout)
                    break
                end
            end

            if all(icheckedout)
                break
            end
            
            j = j + 1;
            polyi{j} = i;
            i_start = i;
        end
    end
end

% Plot for an intermediate check
if plot_level >= 2
    notco = find(icheckedout == 0);
    figure
    triplot(mch.t,mch.p(:,1),mch.p(:,2),'-')
    hold on
    plot(mch.p(polyi{1},1),mch.p(polyi{1},2),'-',LineWidth=2)
    plot(mch.p(notco,1),mch.p(notco,2),'s')
    axis equal
end

tPrev = toc_disp(tPrev,'16-2');

bndpolych = polyi;

% Plot for an intermediate check
if plot_level >= 2
    figure
    for i=1:length(bndpolych)
        poly = bndpolych{i};
        polyxy = mch.coord(poly,:);
        plot(polyxy(:,1),polyxy(:,2),'o-',LineWidth=2)
        hold on
    end
end

% - Mainland bondaries and IBTYPE=64 weir pairs
landbnds = {};
flowbnds = {};
ib64bnds = {};
landbnd = [];
flowbnd = [];
ib64bnd = [];
for ipoly=1:length(bndpolych)
    poly = bndpolych{ipoly};
    for i=1:length(poly)
        ind = poly(i);
        lonlat_ch = mch.p(ind,:);
        jnd = dsearchn(mconn.p,double(mconn.t),lonlat_ch);
        lonlat_cn = mconn.p(jnd,:);
        % Check to see if there is a connecting mesh node which is exactly 
        % at the same location as the channel node. If true, then the
        % node pair is along an IBTYPE=64 boundary. If not, the channel
        % node belongs to a land boundary and/or flow boundary.
        if m_lldist([lonlat_ch(1) lonlat_cn(1)],[lonlat_ch(2) lonlat_cn(2)])*1000.0 < 0.1
            if ~isempty(landbnd)
                landbnd(end+1) = ind;
                landbnds{end+1} = landbnd;
                landbnd = [];
            end
%             if ~isempty(flowbnd)
%                 flowbnds{end+1} = flowbnd;
%                 flowbnd = [];
%             end
            ib64bnd(end+1,1:2) = [ind jnd];
        else
            if ~isempty(ib64bnd)
                ib64bnds{end+1} = flip(ib64bnd);
                landbnd(end+1) = ib64bnd(end,1);
                ib64bnd = [];
            end
            % Check to see if the channel node is one of the end point
            % nodes
            if any(lonlat_ch(1) == mch_lonlat_endpoints1(:,1) & ...
                   lonlat_ch(2) == mch_lonlat_endpoints1(:,2)) || ...
               any(lonlat_ch(1) == mch_lonlat_endpoints2(:,1) & ...
                   lonlat_ch(2) == mch_lonlat_endpoints2(:,2))              
                if ~isempty(landbnd)
                    landbnd(end+1) = ind;
                    landbnds{end+1} = landbnd;
                    landbnd = [];
                end
            else
                landbnd(end+1) = ind;
            end
        end
        % Check to see if the channel node is one of the end point
        % nodes
        if any(lonlat_ch(1) == mch_lonlat_endpoints1(:,1) & ...
               lonlat_ch(2) == mch_lonlat_endpoints1(:,2)) || ...
           any(lonlat_ch(1) == mch_lonlat_endpoints2(:,1) & ...
               lonlat_ch(2) == mch_lonlat_endpoints2(:,2))              
            flowbnd(end+1) = ind;
        else
            if ~isempty(flowbnd)
                flowbnds{end+1} = flowbnd;
                flowbnd = [];
            end
        end
    end
    if ~isempty(landbnd)
        if length(landbnd) >= 2
            landbnds{end+1} = landbnd;
        end
        landbnd = [];
    end
    if ~isempty(flowbnd)
        if length(flowbnd) >= 2
            flowbnds{end+1} = flowbnd;
        end
        flowbnd = [];
    end
    if ~isempty(ib64bnd)
        if length(ib64bnd) >= 2
            ib64bnds{end+1} = flip(ib64bnd);
        end
        ib64bnd = [];
    end
end

% Plot for an intermediate check
if plot_level >= 2
    figure;
    triplot(mconnch.t,mconnch.p(:,1),mconnch.p(:,2))
    hold on
    for i=1:length(landbnds)
        bnd = landbnds{i};
        plot(mch.p(bnd,1),mch.p(bnd,2),'g.-',LineWidth=2)
    end
    for i=1:length(flowbnds)
        bnd = flowbnds{i};
        plot(mch.p(bnd,1),mch.p(bnd,2),'c.-',LineWidth=2)
    end
    for i=1:length(ib64bnds)
        bnd = ib64bnds{i};
        plot(mch.p(bnd(:,1),1),mch.p(bnd(:,1),2),'b.-',LineWidth=2)
        plot(mconn.p(bnd(:,2),1),mconn.p(bnd(:,2),2),'r.-',LineWidth=1)
    end
    axis equal
end

tPrev = toc_disp(tPrev,'17');

if ~no_channelmesh
    % Set the boundaries on mconnch
    
    % - Initialize
    mconnch.bd.nbou   = 0;
    mconnch.bd.nvel   = 0;
    mconnch.bd.nvell  = [];
    mconnch.bd.ibtype = [];
    mconnch.bd.nbvv   = sparse(0,0);
    mconnch.bd.ibconn = sparse(0,0);
    mconnch.bd.barinht= sparse(0,0);
    mconnch.bd.barincfsb = sparse(0,0);
    mconnch.bd.barincfsp = sparse(0,0);
    nnconn = size(mconn.p,1);
    
    % - Set the land boundaries
    for j=1:length(landbnds)
        bnd = landbnds{j};
        nbnd = length(bnd);
        ibd = mconnch.bd.nbou + 1;
        mconnch.bd.nbou = ibd;
        mconnch.bd.nvel = mconnch.bd.nvel + nbnd;
        mconnch.bd.nvell(ibd)  = nbnd;
        mconnch.bd.ibtype(ibd) = 20;
        mconnch.bd.nbvv(1:nbnd,ibd) = bnd+nnconn;    
    end
    
    % - Set the flow boundaries
    for j=1:length(flowbnds)
        bnd = flowbnds{j};
        nbnd = length(bnd);
        ibd = mconnch.bd.nbou + 1;
        mconnch.bd.nbou = ibd;
        mconnch.bd.nvel = mconnch.bd.nvel + nbnd;
        mconnch.bd.nvell(ibd)  = nbnd;
        mconnch.bd.ibtype(ibd) = 22;
        mconnch.bd.nbvv(1:nbnd,ibd) = bnd+nnconn;    
    end
    
    % - Set the ibtype=64 (VEW) boundaries
    for j=1:length(ib64bnds)
        bnd = ib64bnds{j};
        bndch = bnd(:,1); % channel mesh
        bndcn = bnd(:,2); % connecting mesh
        nbnd = size(bnd,1);
        ibd = mconnch.bd.nbou + 1;
        mconnch.bd.nbou = ibd;
        mconnch.bd.nvel = mconnch.bd.nvel + nbnd;
        mconnch.bd.nvell(ibd)  = nbnd;
        mconnch.bd.ibtype(ibd) = 64;
        mconnch.bd.nbvv(1:nbnd,ibd) = bndcn;   
        mconnch.bd.ibconn(1:nbnd,ibd) = bndch+nnconn;
        mconnch.bd.barinht(1:nbnd,ibd)   = -min(mch.b(bndch),mconn.b(bndcn))+1e-4;
        mconnch.bd.barincfsb(1:nbnd,ibd) = 1.0;
        mconnch.bd.barincfsp(1:nbnd,ibd) = 1.0;
    end
    
    % Set coords
    [msub.coord(1:nnsub,1),msub.coord(1:nnsub,2)] ...
        = deg2cart(msub.p(:,2),msub.p(:,1),lat00,lon00);
    [mconnch.coord(1:nnconnch,1),mconnch.coord(1:nnconnch,2)] ...
        = deg2cart(mconnch.p(:,2),mconnch.p(:,1),lat00,lon00);
end

if write_output
    mconnch.write([out14 '_mconnch'],{'13','14'})
end

% Merge the connection+channel mesh to the original mesh

tPrev = toc_disp(tPrev,'18');

% Remove the overlapping elements
% - Update nd_in_removed
nd_in_removed_org = nd_in_removed;
nd_in_removed = union(nd_in_removed_org,setdiff(nd_boundary0,nd_boundary));
% - Find nodes inside polygons and add them
nd_in_polys = [];
nd_on_polys = [];
polynodes = [];
for j=1:length(polys_all)
    poly = polys_all{j};
    polynodes = [polynodes; poly];
    [in,on] = inpolygon(m.coord(:,1),m.coord(:,2),poly(:,1),poly(:,2));
    nd_in_polys = [nd_in_polys;find(in)];
    nd_on_polys = [nd_on_polys;find(on)];
end
nd_in_polys = unique(nd_in_polys);
nd_on_polys = unique(nd_on_polys);
nd_in_polys = setdiff(nd_in_polys, nd_boundary_poly);
nd_in_removed = union(nd_in_removed, nd_in_polys);
% % - Find elements inside polygons and add them
% elems_in_polys = [];
% elems_on_polys = [];
% lons = m.coord(:,1);
% lats = m.coord(:,2);
% elems_lon = mean(lons(m.t(:,:)),2);
% elems_lat = mean(lats(m.t(:,:)),2);
% polynodes = [];
% for j=1:length(polys_all)
%     poly = polys_all{j};
%     polynodes = [polynodes; poly];
%     [in,on] = inpolygon(elems_lon,elems_lat,poly(:,1),poly(:,2));
%     elems_in_polys = [elems_in_polys;find(in)];
%     elems_on_polys = [elems_on_polys;find(on)];
% end
% polynodes = unique(polynodes);
% elems_in_polys = unique(elems_in_polys);
% elems_on_polys = unique(elems_on_polys);
% - Make n2e
nn = size(m.p,1);
ne = size(m.t,1);
n2e = cell(nn,1);
for j = 1:ne
    for nm = 1:3
        n2e{m.t(j,nm)} = [n2e{m.t(j,nm)} j];
    end
end
for i = 1:nn
    n2e{i} = unique(n2e{i});
end
% - Make el_toberemoved
el_toberemoved = [];
for i = nd_in_removed'
    el_toberemoved = [el_toberemoved n2e{i}];
end
% el_toberemoved = union(el_toberemoved, elems_in_polys);

% Find hanging nodes
% - Calculate adjacency matrix of t
% S = sparse(ma.t(:,[1,1,2,2,3,3]),ma.t(:,[2,3,1,3,1,2]),1,nnma,nnma);
% W = sum(S,2);
% hanging_node = find(W==0);
% hanging_nodes = m.t(elems_in_polys,:);
% hanging_nodes = hanging_nodes(:);
% hanging_nodes = setdiff(hanging_nodes, polynodes);
% nd_in_removed = union(nd_in_removed, hanging_nodes);
% Plot for an intermediate check
if plot_level >= 2
    figure;
    triplot(m.t,m.p(:,1),m.p(:,2))
    hold on
    triplot(m.t(el_toberemoved,:),m.p(:,1),m.p(:,2),Color='red')
    plot(m.p(nd_in_removed,1),m.p(nd_in_removed,2),'o',LineWidth=2)
    plot(m.p(setdiff(nd_boundary0,nd_boundary),1),m.p(setdiff(nd_boundary0,nd_boundary),2),'o',LineWidth=1)
    plot(m.p(:,1),m.p(:,2),'k.')
    axis equal
end

% - Remove
ma = m;
nnma = size(ma.p,1);
ma.t(el_toberemoved,:) = [];
if ~isfield(ma.bd,'ibconn')
    ma.bd.ibconn = sparse(size(ma.bd.nbvv,1),size(ma.bd.nbvv,2));
    ma.bd.barinht   = sparse(size(ma.bd.nbvv,1),size(ma.bd.nbvv,2));
    ma.bd.barincfsb = sparse(size(ma.bd.nbvv,1),size(ma.bd.nbvv,2));
    ma.bd.barincfsp = sparse(size(ma.bd.nbvv,1),size(ma.bd.nbvv,2));
end

% Remove the overlapping nodes
mapnd = [(1:nnma)' ones(nnma,1)];
mapnd(nd_in_removed,2) = 0;
mapnd(mapnd(:,2)==0,:) = [];
mapnd(:,2) = 1:size(mapnd,1);  % map from ma to m

mapnd2 = [(1:nnma)' ones(nnma,1)*-1];
mapnd2(mapnd(:,1),2) = mapnd(:,2);  % map from m to ma

ma.p(nd_in_removed,:) = [];
ma.b(nd_in_removed,:) = [];

nnma = size(ma.p,1);

% Update the element table
ma.t(:,1) = mapnd2(ma.t(:,1),2);
ma.t(:,2) = mapnd2(ma.t(:,2),2);
ma.t(:,3) = mapnd2(ma.t(:,3),2);

% Update the boundary node IDs
inonzero = ma.op.nbdv > 0;
ma.op.nbdv(inonzero) = mapnd2(ma.op.nbdv(inonzero),2);
assert(sum(ma.op.nbdv < 0) == 0)

inonzero = ma.bd.nbvv > 0;
ma.bd.nbvv(inonzero) = mapnd2(ma.bd.nbvv(inonzero),2);

inonzero = ma.bd.ibconn > 0;
ma.bd.ibconn(inonzero) = mapnd2(ma.bd.ibconn(inonzero),2);

i = 1;
while true
    doitagain = false;
    for j=1:ma.bd.nvell(i)
        if ma.bd.nbvv(j,i) < 0
            if j > 2 && j < ma.bd.nvell(i)-1
                ibtype2 = ma.bd.ibtype(i);
                nbvv2 = ma.bd.nbvv(j+1:end,i);
                ibconn2 = ma.bd.ibconn(j+1:end,i);
                barinht2 = ma.bd.barinht(j+1:end,i);
                barincfsb2 = ma.bd.barincfsb(j+1:end,i);
                barincfsp2 = ma.bd.barincfsp(j+1:end,i);

                ma.bd.nvell(i) = j-1;
                ma.bd.nbvv(j:end,i) = 0;
                ma.bd.ibconn(j:end,i) = 0;
                ma.bd.barinht(j:end,i) = 0;
                ma.bd.barincfsb(j:end,i) = 0;
                ma.bd.barincfsp(j:end,i) = 0;

                ma.bd.nvell(ma.bd.nbou+1) = length(nonzeros(nbvv2));
                ma.bd.ibtype(ma.bd.nbou+1) = ibtype2;
                ma.bd.nbvv(1:length(nbvv2),ma.bd.nbou+1) = nbvv2;
                ma.bd.ibconn(1:length(nbvv2),ma.bd.nbou+1) = ibconn2;
                ma.bd.barinht(1:length(nbvv2),ma.bd.nbou+1) = barinht2;
                ma.bd.barincfsb(1:length(nbvv2),ma.bd.nbou+1) = barincfsb2;
                ma.bd.barincfsp(1:length(nbvv2),ma.bd.nbou+1) = barincfsp2;

                ma.bd.nbou = ma.bd.nbou + 1;

            elseif j <= 2 && j < ma.bd.nvell(i)-1
                nbvv2 = ma.bd.nbvv(j+1:end,i);
                ibconn2 = ma.bd.ibconn(j+1:end,i);
                barinht2 = ma.bd.barinht(j+1:end,i);
                barincfsb2 = ma.bd.barincfsb(j+1:end,i);
                barincfsp2 = ma.bd.barincfsp(j+1:end,i);

                ma.bd.nbvv(:,i) = 0;
                ma.bd.ibconn(:,i) = 0;
                ma.bd.barinht(:,i) = 0;
                ma.bd.barincfsb(:,i) = 0;
                ma.bd.barincfsp(:,i) = 0;

                ma.bd.nvell(i) = length(nonzeros(nbvv2));
                ma.bd.nbvv(1:length(nbvv2),i) = nbvv2;
                ma.bd.ibconn(1:length(nbvv2),i) = ibconn2;
                ma.bd.barinht(1:length(nbvv2),i) = barinht2;
                ma.bd.barincfsb(1:length(nbvv2),i) = barincfsb2;
                ma.bd.barincfsp(1:length(nbvv2),i) = barincfsp2;

                doitagain = true;

            elseif j > 2 && j >= ma.bd.nvell(i)-1
                ma.bd.nvell(i) = j-1;
                ma.bd.nbvv(j:end,i) = 0;
                ma.bd.ibconn(j:end,i) = 0;
                ma.bd.barinht(j:end,i) = 0;
                ma.bd.barincfsb(j:end,i) = 0;
                ma.bd.barincfsp(j:end,i) = 0;

            else
                ma.bd.nvell(i) = [];
                ma.bd.ibtype(i) = [];
                ma.bd.nbvv(:,i) = [];
                ma.bd.ibconn(:,i) = [];
                ma.bd.barinht(:,i) = [];
                ma.bd.barincfsb(:,i) = [];
                ma.bd.barincfsp(:,i) = [];

                ma.bd.nbou = ma.bd.nbou - 1;

                doitagain = true;
            end
            break
        end
    end
    if ~doitagain
        if i >= ma.bd.nbou
            break
        end
        i = i + 1;
    end
end
ma.bd.nvel = 0;
for i=1:ma.bd.nbou
    ma.bd.nvel = ma.bd.nvel + ma.bd.nvell(i);
end

% Plot for an intermediate check
if plot_level >= 2
    figure
    triplot(ma.t,ma.p(:,1),ma.p(:,2))
    axis equal
end

% Remove nodal attributes at the removed nodes and remap

for k=1:ma.f13.nAttr
    nds = ma.f13.userval.Atr(k).Val(1,:);
    ir = ismember(nds,nd_in_removed);
    ma.f13.userval.Atr(k).Val(:,ir) = [];
    ma.f13.userval.Atr(k).Val(1,:) = ...
        mapnd2(ma.f13.userval.Atr(k).Val(1,:),2);
    ma.f13.userval.Atr(k).usernumnodes ...
        = size(ma.f13.userval.Atr(k).Val,2);
end

% Find the corresponding nodes of ma to the nodes in mconnch
mconnch2ma = zeros(nnconnch,2); % map from mconnch to ma
mconnch2ma(:,1) = 1:nnconnch;
for i=1:nnconnch
    px = mconnch.p(i,1);
    py = mconnch.p(i,2);
    ima = find(ma.p(:,1) == px & ma.p(:,2)==py);
    if ~isempty(ima)
        mconnch2ma(i,2) = ima;
    end
end
nnconnch_stay = sum(mconnch2ma(:,2) == 0);
mconnch2ma(mconnch2ma(:,2) == 0,2) = (1:nnconnch_stay)'+nnma;
map = mconnch2ma;

% map = [(1:nnconnch)' (1:nnconnch)'+nnma-(nn_outer-1)];
% for i = 1:nn_outer-1
%     p.x = mconnch.p(i,1);
%     p.y = mconnch.p(i,2);
%     ima = find(ma.p(:,1) == p.x & ma.p(:,2)==p.y);
%     map(i,2) = ima;
% end

tPrev = toc_disp(tPrev,'19');

% map from msub to ma
map_msub2ma = mapnd2(nd_in,2);

% Shift for the node IDs of the conn+ch mesh
mb = mconnch;
mb.t(:,1) = map(mb.t(:,1),2);
mb.t(:,2) = map(mb.t(:,2),2);
mb.t(:,3) = map(mb.t(:,3),2);

mb.p(map(:,2)<=nnma,:) = [];
mb.b(map(:,2)<=nnma) = [];

if isfield(mb.bd,'nbvv')
    nonzero = mb.bd.nbvv ~= 0;
    mb.bd.nbvv(nonzero) = map(mb.bd.nbvv(nonzero),2);
end
if isfield(mb.bd,'ibconn')
    nonzero = mb.bd.ibconn ~= 0;
    mb.bd.ibconn(nonzero) = map(mb.bd.ibconn(nonzero),2);
end

[xs,ys] ...
    = deg2cart(mb.p(:,2),mb.p(:,1),lat00,lon00);
mb.coord = [xs ys];

ma_org = ma;

% Add mb to ma
ma.p = [ma.p;mb.p];
ma.b = [ma.b;mb.b];
ma.t = [ma.t;mb.t];

if ~no_channelmesh
    nbou0 = ma.bd.nbou;
    ma.bd.nbou = ma.bd.nbou + mb.bd.nbou;
    ma.bd.nvel = ma.bd.nvel + mb.bd.nvel;
    ma.bd.nvell = [ma.bd.nvell mb.bd.nvell];
    ma.bd.ibtype = [ma.bd.ibtype mb.bd.ibtype];
    ma.bd.nbvv(1:size(mb.bd.nbvv,1),(nbou0+1):(nbou0+mb.bd.nbou)) ...
        = mb.bd.nbvv;
    ma.bd.ibconn(1:size(mb.bd.ibconn,1),(nbou0+1):(nbou0+mb.bd.nbou)) ...
        = mb.bd.ibconn;
    ma.bd.barinht(1:size(mb.bd.ibconn,1),(nbou0+1):(nbou0+mb.bd.nbou)) ...
        = mb.bd.barinht;
    ma.bd.barincfsb(1:size(mb.bd.ibconn,1),(nbou0+1):(nbou0+mb.bd.nbou)) ...
        = mb.bd.barincfsb;
    ma.bd.barincfsp(1:size(mb.bd.ibconn,1),(nbou0+1):(nbou0+mb.bd.nbou)) ...
        = mb.bd.barincfsp;
end

% Inerpolate nodal attribute values at nodes of mb part
mapnd22 = mapnd2(:,2);
mapnd22(end+1) = -1;
map2 = map(:,2);
map2(end+1) = -1;
mapb = map(map(:,2)>nnma,:); % mapb = map(map(:,2)>=nnma,:);
nnmb = size(mb.p,1);
tri = msub.t;
N = [];
tinan = [];
for k=1:msub.f13.nAttr
    Atr = msub.f13.defval.Atr(k);
    defvals = Atr.Val';
    nvals = Atr.ValuesPerNode;
    
    Atr = msub.f13.userval.Atr(k);
    uservals = Atr.Val';

    if strcmp(Atr.AttrName,'condensed_nodes')
        val1 = ma.f13.userval.Atr(k).Val(1,:); % already mapped to ma

        val2 = ma.f13.userval.Atr(k).Val(2:end,:);
        val2(val2==0) = length(mapnd22);
        val2 = mapnd22(val2); % map from m to ma
        val2(val2==-1) = 0;
        if size(val1,2) ~= size(val2,2)
            val2 = val2';
        end

        val12 = [val1;val2];

        val3 = msub.f13.userval.Atr(k).Val(1,:);
        val3 = map_msub2ma(val3); % map msub to ma

        val4 = msub.f13.userval.Atr(k).Val(2:end,:);
        val4(val4==0) = length(mapnd22); % map m to ma
        val4 = mapnd22(val4);
        val4(val4==-1) = 0;
        if size(val3,1) ~= 1
            val3 = val3';
        end
        if size(val3,2) ~= size(val4,2)
            val4 = val4';
        end

        val34 = [val3;val4];
        Val = [val12 val34];

%         val3 = msub.f13.userval.Atr(k).Val;
%         val3(val3==0) = length(map2);
%         val3 = map2(val3);
%         val3(val3==-1) = 0;
%         Val = [val12 val3];
        ma.f13.userval.Atr(k).Val = Val;
        ma.f13.userval.Atr(k).usernumnodes = size(ma.f13.userval.Atr(k).Val,2);
    else
        fullvals = repmat(defvals,nnsub,1);
        fullvals(uservals(:,1),:) = uservals(:,2:nvals+1);
        
        % fullvals = make_vals_for_interp(msub_nodes_inchannel,fullvals,neinodessub);
        
        [zi,N,tri,tinan] = interptri_multi(tri,msub.coord(:,1),msub.coord(:,2),fullvals,mb.coord(:,1),mb.coord(:,2),N,tinan);
        
        idxnan = isnan(zi);
        zi(idxnan) = defvals(1);
        
        idxdef = sum(zi == defvals,2) == nvals;
        zi = [mapb(:,2) zi];
        zi(idxdef,:) = [];
        nuservals2 = size(zi,1);
        
        ma.f13.userval.Atr(k).Val(:,end+1:end+nuservals2) = zi';
        ma.f13.userval.Atr(k).usernumnodes = size(ma.f13.userval.Atr(k).Val,2);
    end
end
ma.f13.NumOfNodes = size(ma.p,1);

tPrev = toc_disp(tPrev,'20');

flags = zeros(size(ma.p,1),1);
flags(ma.t(:,1)) = 1;
flags(ma.t(:,2)) = 1;
flags(ma.t(:,3)) = 1;
inds = find(flags==0);

% Plot for an intermediate check
if plot_level >= 2
    figure
    triplot(ma.t, ma.p(:,1), ma.p(:,2), 'k-')
    hold on
    plot(ma.p(inds,1), ma.p(inds,2), 'o')
    axis equal
    
    figure
    triplot(mconn.t, mconn.p(:,1), mconn.p(:,2), 'k-')
    hold on
    plot(ma.p(inds,1), ma.p(inds,2), 'o')
    axis equal
    
    figure
    triplot(ma_org.t, ma_org.p(:,1), ma_org.p(:,2), 'k-')
    hold on
    plot(ma.p(inds,1), ma.p(inds,2), 'o')
    axis equal
    
    figure
    triplot(msub.t, msub.p(:,1), msub.p(:,2), 'k-')
    hold on
    plot(ma.p(inds,1), ma.p(inds,2), 'o')
    axis equal
end

% Renum
if do_renum
    [ma, map_renum] = renum_mesh_and_condensednodes(ma);
else
    map_renum = 1:length(ma.b);
end

if ~no_channelmesh
    % - map channel bank nodes
    for i=1:shp.nseq
        shp.seq(i).bn.gid ...
            = map_node(shp.seq(i).bn.gid);
    end
end

if ~no_channelmesh
    % Build nodal attributes including condensed_nodes
    % - CONDENSED_NODES
    condensed_nodes = [];
    channel_nodes = [];
    
    % -- Gather cross point ids
    for i=1:length(shp.endpoints.connected_points)
        crossp = shp.endpoints.connected_points(i).crossp;
        ncrossp = length(crossp);
        if ncrossp == 2
            if crossp(1).removed == 0 && crossp(2).removed == 0
                ncn = length(crossp(2).mid_gid) + 2;
                channel_nodes(end+1:end+ncn) = [ ...
                    map_renum(map2(map_remove_collapsed2(crossp(1).gid)+nnconn)), ...
                    map_renum(map2(map_remove_collapsed2(crossp(2).gid)+nnconn)), ...
                    map_renum(map2(map_remove_collapsed2(crossp(2).mid_gid')+nnconn)), ...
                    ];
                if isempty(crossp(2).mid_gid)
                    condensed_nodes(end+1,1) = 0;
                    condensed_nodes(end,1) = map_renum(map2(map_remove_collapsed2(crossp(1).gid)+nnconn));
                    condensed_nodes(end,2) = map_renum(map2(map_remove_collapsed2(crossp(2).gid)+nnconn));
                end
            end
        elseif ncrossp >=3
            ncn = length(shp.endpoints.connected_points(i).mesh.gid);
            channel_nodes(end+1:end+ncn) = ...
                map_renum(map2(map_remove_collapsed2( ...
                   shp.endpoints.connected_points(i).mesh.gid)+nnconn ...
                ));

            min_esize = 1e12;
            for j=1:size(shp.endpoints.connected_points(i).mesh.elem,1)
                el = shp.endpoints.connected_points(i).mesh.elem(j,:);
                xys = [shp.endpoints.connected_points(i).mesh.coord(el,1), ...
                       shp.endpoints.connected_points(i).mesh.coord(el,2)];
                min_esize = min(min_esize, calculate_element_size(xys));
            end
            if min_esize < min_spacing*0.5 || ...
               (min_esize < min_spacing && ncn == 3)
                condensed_nodes(end+1,1) = 0;
                % for j=1:ncrossp
                %     condensed_nodes(end,j) = map_renum(map2(map_remove_collapsed2(crossp(j).gid)+nnconn));
                % end
                condensed_nodes(end,1:ncn) = ...
                    map_renum(map2(map_remove_collapsed2( ...
                       shp.endpoints.connected_points(i).mesh.gid) + nnconn ...
                    ));
            end
        end
    end

    % -- Gather channel node ids
    gid_crossp = [];
    for i=1:length(shp.endpoints.connected_points)
        for j=1:length(shp.endpoints.connected_points(i).crossp)
            gid_crossp(end+1) = shp.endpoints.connected_points(i).crossp(j).gid;
        end
    end
    for i=1:shp.nseq
        % n = length(shp.seq(i).mesh.gid);
        % ista = 3;
        % iend = n - 3;
        % condensed_nodes(end+1:end+n/2-2,1:2) = map_renum(map2([map_remove_collapsed2(shp.seq(i).mesh.gid(ista:2:iend))+nnconn map_remove_collapsed2(shp.seq(i).mesh.gid(ista+1:2:iend+1))+nnconn]));
        mesh_gid = shp.seq(i).mesh.gid;
        mesh_gid(mesh_gid == 0) = [];
        ncn = length(mesh_gid);
        channel_nodes(end+1:end+ncn) = ...
            map_renum(map2(map_remove_collapsed2(mesh_gid)+nnconn));
        
        not_splitted = find(~shp.seq(i).cl.splitted);
        % not_splitted(ismember(not_splitted,[1,size(shp.seq(i).cl.coord,1)])) = [];
        if not_splitted(1) == 1 && ismember(shp.seq(i).bn.gid(1), gid_crossp)
            not_splitted(1) = [];
        end
        if not_splitted(end) == length(~shp.seq(i).cl.splitted) && ismember(shp.seq(i).bn.gid(end), gid_crossp)
            not_splitted(end) = [];
        end
        bnds = [(not_splitted*2-1)', (not_splitted*2)'];
        xy1 = [shp.seq(i).bn.coord(bnds(:,1),1), shp.seq(i).bn.coord(bnds(:,1),2)];
        xy2 = [shp.seq(i).bn.coord(bnds(:,2),1), shp.seq(i).bn.coord(bnds(:,2),2)];
        dxy = xy1 - xy2;
        len = sqrt(dxy(:,1).^2 +dxy(:,2).^2);
        bnds(len > pave_opts.min_condensed_nodes_spacing,:) = [];
        bnd_gids = [shp.seq(i).bn.gid(bnds(:,1)) shp.seq(i).bn.gid(bnds(:,2))];
        if ~isempty(bnd_gids)
            bnd_gids(bnd_gids(:,1) == 0 | bnd_gids(:,2) == 0, :) = []; % drop collapsed nodes
            condensed_nodes(end+1:end+size(bnd_gids,1),1:2) = map_renum(map2([map_remove_collapsed2(bnd_gids(:,1))+nnconn map_remove_collapsed2(bnd_gids(:,2))+nnconn]));
        end
    end

    % -- Add condensed_nodes
    ma = add_condensed_nodes_to_f13(ma,condensed_nodes);

    % - Set channel values
    if ~isempty(default_ndattrs) && ~isempty(channel_ndattrs)
        ma = set_channel_nodal_attributes(ma,channel_nodes,default_ndattrs,channel_ndattrs);
    end
end

% Plot for an intermediate check
if plot_level >= 2
    figure
    triplot(ma.t, ma.p(:,1), ma.p(:,2))
    hold on
    plot(ma.p(channel_nodes,1), ma.p(channel_nodes,2), 'ro')
    for i=channel_nodes
        text(ma.p(i,1),ma.p(i,2),int2str(i))
    end
    axis equal
end

% Sort the nodal attribute values according to the node id
for k=1:ma.f13.nAttr
    Atr = ma.f13.defval.Atr(k);
    nvals = Atr.ValuesPerNode;
    Val = ma.f13.userval.Atr(k).Val;
    [sorted,I] = sort(Val(1,:));
    Val = Val(:,I);
    ma.f13.userval.Atr(k).Val = Val;
end

nnma = size(ma.p,1);
nema = size(ma.t,1);
% if ~no_channelmesh
%     % Make land boundaries for the remaining parts of the domain boundary
%     % - Make n2e
%     n2ema = cell(nnma,1);
%     for j = 1:nema
%         for nm = 1:3
%             n2ema{ma.t(j,nm)} = [n2ema{ma.t(j,nm)} j];
%         end
%     end
%     for i = 1:nnma
%         n2ema{i} = unique(n2ema{i});
%     end
%     % - Make neinodes
%     neinodesma = cell(nnma,1);
%     for i = 1:nnma
%         neielems = n2ema{i,1};
%         nei = ma.t(neielems,:);
%         neinodesma{i} = unique(nei(:)');
%     end
%     % - Make a list of nodes along the domain boundary
%     boundary_nodes = [];
%     for i=1:nnma
%         neielems = n2ema{i};
%         neinodes = neinodesma{i};
%         neinodes(neinodes == i) = [];
% 
%         count = zeros(length(neinodes),1);
%         for j=neielems
%             nm = ma.t(j,:);
%             i0 = find(nm == i);
%             i1 = mod(i0-1+1,3)+1;
%             i2 = mod(i0-1+2,3)+1;
%             ii1 = find(neinodes == nm(i1));
%             assert(~isempty(ii1));
%             ii2 = find(neinodes == nm(i2));
%             assert(~isempty(ii2));
%             count(ii1) = count(ii1) + 1;
%             count(ii2) = count(ii2) + 1;
%         end
% 
%         if any(count == 1)
%             boundary_nodes(end+1) = i;
%         end
%     end
%     tPrev = toc_disp(tPrev,'20-1');
%     % - Make a list of the boundary nodes where the boundary conditions are
%     % assigned.
%     bcnodes = [nonzeros(ma.op.nbdv);nonzeros(ma.bd.nbvv);flip(nonzeros(ma.bd.ibconn))];
%     % - Make boundary polylines
%     ibnd_start = 1;
%     ibnd = ibnd_start;
%     j = 1;
%     polyi = {};
%     polyi{j} = boundary_nodes(ibnd); % first node
%     icheckedout = zeros(length(boundary_nodes),1);
%     icheckedout(ibnd) = 1;
%     cnt = 0;
%     seqseq = [boundary_nodes(ibnd)];
%     while true
%         i = boundary_nodes(ibnd);
%         neiele = n2ema{i};
%         found = false;
%         for e = neiele
%             nm = ma.t(e,:);
%             ii = find(nm == i);
%             ii_next = mod(ii+1-1,3)+1;
%             i_next = nm(ii_next);
%             neiele_next = n2ema{i_next};
%             found2 = false;
%             for ee = neiele
%                 if e == ee, continue, end
%                 nm = ma.t(ee,:);
%                 ij = find(nm == i);
%                 ij_next = mod(ij-1-1,3)+1;
%                 j_next = nm(ij_next);
%                 if j_next == i_next
%                     found2 = true;
%                     break
%                 end
%             end
%             if ~found2
%                 assert(ismember(i_next,boundary_nodes))
%                 found = true;
%                 break;
%             end
%         end
% 
%         assert(found)
% 
%         polyi{j} = [polyi{j} i_next];
%         seqseq = [seqseq, i_next];
%         ibnd_next = find(boundary_nodes == i_next);
%         assert(~isempty(ibnd_next))
%         icheckedout(ibnd_next) = 1;
%         ibnd = ibnd_next;
% 
%         if ibnd_next == ibnd_start
%             if sum(icheckedout) == length(icheckedout)
%                 break;
%             else
%                 j = j + 1;
%                 ibnd = find(icheckedout == 0, 1);
%                 polyi{j} = boundary_nodes(ibnd);
%                 icheckedout(ibnd) = 1;
%                 ibnd_start = ibnd;
%                 seqseq = [seqseq, polyi{j}];
%             end
%         end
%         cnt = cnt + 1;
%         if cnt >= 30609*3
%             break;
%         end
%     end
%     bndpolyma = polyi;
%     tPrev = toc_disp(tPrev,'20-2');
%     % figure
%     % hold on
%     % notco = [];
%     % nco = zeros(size(ma.p,1),1);
%     % for i=1:length(polyi)
%     %     polyii = polyi{i};
%     %     plot(ma.p(polyii,1),ma.p(polyii,2),'-')
%     %     notco = [notco, polyii];
%     %     nco(polyii) = nco(polyii) + 1;
%     % end
%     % notco = setdiff(boundary_nodes,notco);
%     % plot(ma.p(notco,1),ma.p(notco,2),'s')
%     % plot(ma.p(nco>1,1),ma.p(nco>1,2),'o')
%     % plot(ma.p(seqseq(end-10:end),1),ma.p(seqseq(end-10:end),2),'-',LineWidth=2)
%     % plot(ma.p(seqseq(end),1),ma.p(seqseq(end),2),'o',LineWidth=2)
%     % axis equal
%     % - Assign land boundaries along the sequences of nodes where no boundary
%     % condition is assigned.
%     landbnds = {};
%     landbnd = [];
%     for i=1:length(bndpolyma)
%         poly = bndpolyma{i};
%         for j=1:length(poly)
%             idx = find(bcnodes==poly(j));
%             if poly(j) == 5117
%                 disp([i j])
%             end
%             onbc = false;
%             if ~isempty(idx) && length(idx) == 1
%                 if j > 1 && idx > 1 && ~isempty(landbnd)
%                     if poly(j-1) == bcnodes(idx-1)
%                         onbc = true;
%                     end
%                 end
%                 if j < length(poly) && idx < length(bcnodes) && isempty(landbnd)
%                     if poly(j+1) == bcnodes(idx+1)
%                         onbc = true;
%                     end
%                 end
%             elseif length(idx) >= 2
%                 onbc = true;
%             end
%             if isempty(idx) || ~onbc
%                 if j > 1 && isempty(landbnd) && isempty(idx)
%                     landbnd = [poly(j-1)];
%                 end
%                 landbnd(end+1) = poly(j);
%             else
%                 if ~isempty(landbnd)
%                     if isempty(idx)
%                         landbnd(end+1) = poly(j);
%                     end
%                     if length(landbnd) >= 2
%                         landbnds{end+1} = landbnd;
%                     end
%                     landbnd = [];
%                 end
%             end
%         end
%         if length(landbnd) >= 2
%             landbnds{end+1} = landbnd;
%         end
%         landbnd = [];
%     end
%     % Plot for an intermediate check
%     if plot_level >= 2
%         figure
%         triplot(ma.t,ma.p(:,1),ma.p(:,2))
%         hold on
%         for i=1:length(landbnds)
%             landbnd = landbnds{i};
%             plot(ma.p(landbnd,1),ma.p(landbnd,2),'o-',LineWidth=2)
%         end
%         plot(ma.p(bcnodes,1),ma.p(bcnodes,2),'ko',LineWidth=2)
%         axis equal
%     end
%     % - Set land BCs
%     for i=1:length(landbnds)
%         landbnd = landbnds{i};
%         nbnd = length(landbnd);
%         ma.bd.nbou = ma.bd.nbou + 1;
%         ma.bd.nvel = ma.bd.nvel + nbnd;
%         ma.bd.nvell(end+1) = nbnd;
%         ma.bd.ibtype(end+1) = 20; % mainland boundary, natural
%         ma.bd.nbvv(1:nbnd,ma.bd.nbou) = landbnd;
%         ma.bd.ibconn(1,ma.bd.nbou) = 0;
%         ma.bd.barinht(1,ma.bd.nbou) = 0;
%         ma.bd.barincfsb(1,ma.bd.nbou) = 0;
%         ma.bd.barincfsp(1,ma.bd.nbou) = 0;
%     end
% end

% Adjust the depths along channels
for i=1:ma.bd.nbou
    if ma.bd.ibtype(i) ~= 64
        continue
    end
    % The depth along a channel should be lower than the pairing floodplain node
    ma.b(nonzeros(ma.bd.ibconn(:,i))) ...
        = max(ma.b(nonzeros(ma.bd.ibconn(:,i))), ...
              ma.b(nonzeros(ma.bd.nbvv(:,i)))+pave_opts.min_channel_depth);
end
for i=1:ma.f13.nAttr
    if strcmp(ma.f13.userval.Atr(i).AttrName,'condensed_nodes')
        for j=1:ma.f13.userval.Atr(i).usernumnodes
            idx = nonzeros(ma.f13.userval.Atr(i).Val(:,j));
            ma.b(idx) = max(ma.b(idx));
        end
    end
end
if ~no_channelmesh
    ma.b(channel_nodes) = min(ma.b(channel_nodes), pave_opts.max_channel_depth);
    for i=1:ma.bd.nbou
        if ma.bd.ibtype(i) ~= 64
            continue
        end
        ma.b(nonzeros(ma.bd.ibconn(:,i))) ...
            = min(ma.b(nonzeros(ma.bd.ibconn(:,i))), ...
                  pave_opts.max_channel_depth);
    end
end

% INITIAL_RIVER_ELEVATION
if ~no_channelmesh && set_initial_river_elevation
    ma = add_initial_river_elevation_to_f13(ma,channel_nodes);
end

% Copy the mesh to the final mesh varialbe
m_updated_global = ma;

tPrev = toc_disp(tPrev,'21');

if plot_level >= 2
    figure
    triplot(m_updated_global.t,m_updated_global.p(:,1),m_updated_global.p(:,2))
    hold on
    axis equal
end


if write_output
    % Write fort.14
    m_updated_global.write(out14,{'13','14'});
    % Rename fort.14
    movefile([out14 '.14'],[out14 '.grd']);
end

tPrev = toc_disp(tPrev,'22');

% % Plot
% m_updated_global.plot('proj','none');

end

function [ma, map_renum] = renum_mesh_and_condensednodes(ma)
[ma,map_renum] = ma.renum();

% Renum condensed_nodes values
for k=1:ma.f13.nAttr
    if strcmp(ma.f13.userval.Atr(k).AttrName,'condensed_nodes')
        val = ma.f13.userval.Atr(k).Val(2:end,:);
        idx = find(val~=0);
        val(idx) = map_renum(val(idx));
        ma.f13.userval.Atr(k).Val(2:end,:) = val;
    end
end
end