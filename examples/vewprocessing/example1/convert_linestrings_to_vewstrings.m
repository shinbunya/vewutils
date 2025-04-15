% Paths to the input/output files
basedir = 'H:\mnt\sshfs\hatteras\home\sbunya\GitHub\adcircutils\adcircutils\vewchannel\examples\example1';
f14file = fullfile(basedir, 'fort.14');
vewlocation_geofile = fullfile(basedir, 'example1_vew_locations.shp');
vewfile = fullfile(basedir, 'vewstings.yaml');

dist_max = 10.0; % Maximum distance for nearest neighbor search in meters

% Read the shapefile
vewlocation_gdf = shaperead(vewlocation_geofile);
disp(vewlocation_gdf);

% Load the mesh file
mesh = msh(f14file); % Custom function to load ADCIRC mesh

x = mesh.nodes(:,1);
y = mesh.nodes(:,2);
tri = mesh.;

% Create KDTree for nearest neighbor search
tree = KDTreeSearcher([x, y]);

nodestrings = {};

for i = 1:length(vewlocation_gdf)
    line = vewlocation_gdf(i).X;
    slx = line(1);
    sly = line(2);
    elx = line(end-1);
    ely = line(end);
    
    [nearest_node_index, distance] = knnsearch(tree, [slx, sly]);
    
    if distance > dist_max
        continue;
    end
    
    nodestring = nearest_node_index;
    
    ipos = 1;
    while ipos < length(line) - 1
        xl0 = line(ipos);
        yl0 = line(ipos+1);
        xl1 = line(ipos+2);
        yl1 = line(ipos+3);
        
        % Find neighboring nodes
        neighbors = find_neighbors(mesh, nearest_node_index); % Custom function
        
        xnei = x(neighbors);
        ynei = y(neighbors);
        
        distances = arrayfun(@(i) distance_point_line([xnei(i), ynei(i)], [xl0, yl0], [xl1, yl1]), 1:length(neighbors));
        [min_distance, min_idx] = min(distances);
        
        if min_distance > dist_max
            ipos = ipos + 2;
            continue;
        end
        
        nearest_node_index = neighbors(min_idx);
        nodestring = [nodestring, nearest_node_index];
        
        if norm([elx, ely] - [x(nearest_node_index), y(nearest_node_index)]) < dist_max
            break;
        end
    end
    
    if norm([elx, ely] - [x(nodestring(1)), y(nodestring(1))]) < dist_max
        nodestring = [nodestring, nodestring(1)];
    end
    
    nodestrings{end+1} = nodestring;
end

% Visualization
figure; hold on;
for i = 1:length(nodestrings)
    plot(x(nodestrings{i}), y(nodestrings{i}), '-o');
end
xlabel('X Coordinate'); ylabel('Y Coordinate');

% Assign values to each node
nodestrings_with_values = cell(size(nodestrings));
for i = 1:length(nodestrings)
    nodestrings_with_values{i} = arrayfun(@(node) struct('node', node, 'x', x(node), 'y', y(node), ...
        'bank_elevation', 1.5 * (y(node) > 500) + -4.0 * (y(node) <= 500), 'bank_mannings_n', 0.03), ...
        nodestrings{i}, 'UniformOutput', false);
end

data = struct('vewstrings', {nodestrings_with_values});
yaml_write(data, vewfile); % Custom function to write YAML