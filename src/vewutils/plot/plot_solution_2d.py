#!/usr/bin/env python3
"""
Plot CG ADCIRC water levels and velocity fields from NetCDF files.
Similar functionality to plot_dg_output.py but for CG ADCIRC output files.
"""

import os
import sys
import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import xarray as xr
from datetime import datetime
import pandas as pd
from scipy.interpolate import griddata
from dateutil import parser as date_parser
import zipfile
import tempfile
import xml.etree.ElementTree as ET

# Cartopy imports for background imagery
try:
    import cartopy.crs as ccrs
    from cartopy.io.img_tiles import OSM
    from cartopy.io import img_tiles as cimgt
    from cartopy import feature as cfeature
    CARTOPY_AVAILABLE = True
except ImportError:
    CARTOPY_AVAILABLE = False
    print("Warning: Cartopy not available. Background imagery will be disabled.")

# Scale bar import
try:
    from matplotlib_scalebar.scalebar import ScaleBar
    SCALEBAR_AVAILABLE = True
except ImportError:
    SCALEBAR_AVAILABLE = False


def to_datetime(date):
    """Convert numpy datetime64 to Python datetime."""
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    return datetime.utcfromtimestamp(timestamp)


def format_datetime_short(dt_str):
    """
    Format datetime string from 'HHMM UTC MMM DD' to 'M/D HHZ' format.
    
    Parameters
    ----------
    dt_str : str
        Datetime string in format like '1200 UTC SEP 14'
        
    Returns
    -------
    str
        Formatted datetime string in 'M/D HHZ' format (e.g., '9/14 12Z'), 
        or original string if parsing fails
    """
    if not dt_str or not dt_str.strip():
        return dt_str
    
    dt_str = dt_str.strip()
    
    # Try to match format: "HHMM UTC MMM DD" (e.g., "1200 UTC SEP 14")
    # Pattern: 4 digits, space, UTC, space, 3-letter month, space, 1-2 digits
    pattern = r'(\d{4})\s+UTC\s+(\w{3})\s+(\d{1,2})'
    match = re.match(pattern, dt_str)
    
    if match:
        hhmm = match.group(1)
        month_str = match.group(2).upper()
        day = int(match.group(3))
        
        # Extract hour from HHMM
        hour = int(hhmm[:2])
        
        # Map month abbreviations to numbers
        month_map = {
            'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
            'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12
        }
        
        if month_str in month_map:
            month = month_map[month_str]
            return f"{month}/{day} {hour:02d}Z"
    
    # If parsing fails, try dateutil parser as fallback
    try:
        dt = date_parser.parse(dt_str)
        return f"{dt.month}/{dt.day} {dt.hour:02d}Z"
    except (ValueError, AttributeError, ImportError):
        pass
    
    # If all parsing fails, return original string
    return dt_str


def get_category_suffix(category_str):
    """
    Convert category string to suffix format for annotation.
    
    Parameters
    ----------
    category_str : str
        Category string (e.g., 'cat1', 'cat2', 'ts', 'ex')
        
    Returns
    -------
    str
        Suffix string: ', TS', ', EX', or ', Cat 1' through ', Cat 5'
    """
    if not category_str:
        return ''
    
    category_lower = category_str.lower().strip()
    
    # Check for storm categories (ts or ex)
    if category_lower in ['ts', 'ex']:
        return f', {category_str.upper()}'
    
    # Check for hurricane categories (cat1, cat2, etc.)
    if category_lower.startswith('cat'):
        try:
            cat_num = int(category_lower[3:])
            if 1 <= cat_num <= 5:
                return f', Cat {cat_num}'
        except ValueError:
            pass
    
    # Try to extract number from category string
    match = re.search(r'cat(\d)', category_lower)
    if match:
        cat_num = int(match.group(1))
        if 1 <= cat_num <= 5:
            return f', Cat {cat_num}'
    
    # If no match, check if it contains 'ts' or 'ex'
    if 'ts' in category_lower:
        return ', TS'
    if 'ex' in category_lower:
        return ', EX'
    
    return ''


def parse_kmz_track(kmz_file):
    """
    Parse a KMZ file and extract hurricane track coordinates, datetimes, and categories.
    
    Parameters
    ----------
    kmz_file : str
        Path to the KMZ file
        
    Returns
    -------
    tuple
        (longitudes, latitudes, datetimes, categories) arrays of track data
    """
    longitudes = []
    latitudes = []
    datetimes = []
    categories = []
    
    try:
        # Extract KMZ file (it's a ZIP file)
        with zipfile.ZipFile(kmz_file, 'r') as kmz:
            # Find the KML file inside
            kml_files = [f for f in kmz.namelist() if f.endswith('.kml')]
            if not kml_files:
                print(f"Warning: No KML file found in {kmz_file}")
                return np.array([]), np.array([]), [], []
            
            # Read the first KML file
            kml_content = kmz.read(kml_files[0])
            
            # Parse XML
            root = ET.fromstring(kml_content)
            
            # Register namespaces
            namespaces = {'kml': 'http://earth.google.com/kml/2.2'}
            
            # Find all Placemark elements
            placemarks = root.findall('.//kml:Placemark', namespaces)
            
            for placemark in placemarks:
                # Extract name (datetime)
                name_elem = placemark.find('.//kml:name', namespaces)
                dt_str = name_elem.text.strip() if name_elem is not None and name_elem.text else ""
                
                # Extract category from styleUrl (e.g., #cat1, #cat2, #ts, #ex)
                category_str = ""
                style_url_elem = placemark.find('.//kml:styleUrl', namespaces)
                if style_url_elem is not None and style_url_elem.text:
                    style_url = style_url_elem.text.strip()
                    # Remove '#' prefix if present
                    if style_url.startswith('#'):
                        category_str = style_url[1:]
                    else:
                        category_str = style_url
                
                # Also check ExtendedData for category if styleUrl doesn't have it
                if not category_str:
                    ext_data = placemark.find('.//kml:ExtendedData', namespaces)
                    if ext_data is not None:
                        for data in ext_data.findall('.//kml:Data', namespaces):
                            name_attr = data.get('name')
                            if name_attr and ('cat' in name_attr.lower() or name_attr.lower() in ['ts', 'ex']):
                                value_elem = data.find('.//kml:value', namespaces)
                                if value_elem is not None and value_elem.text:
                                    category_str = value_elem.text.strip()
                                    break
                
                # Find coordinates within Point elements
                point = placemark.find('.//kml:Point', namespaces)
                if point is not None:
                    coords_elem = point.find('.//kml:coordinates', namespaces)
                    if coords_elem is not None and coords_elem.text:
                        # Parse coordinates: format is "lon,lat,alt"
                        coords = coords_elem.text.strip().split(',')
                        if len(coords) >= 2:
                            try:
                                lon = float(coords[0])
                                lat = float(coords[1])
                                longitudes.append(lon)
                                latitudes.append(lat)
                                datetimes.append(dt_str)
                                categories.append(category_str)
                            except ValueError:
                                continue
                
                # Also check for LineString coordinates (if present)
                linestring = placemark.find('.//kml:LineString', namespaces)
                if linestring is not None:
                    coords_elem = linestring.find('.//kml:coordinates', namespaces)
                    if coords_elem is not None and coords_elem.text:
                        # Parse multiple coordinates (space-separated)
                        coord_strings = coords_elem.text.strip().split()
                        for coord_str in coord_strings:
                            coords = coord_str.split(',')
                            if len(coords) >= 2:
                                try:
                                    lon = float(coords[0])
                                    lat = float(coords[1])
                                    longitudes.append(lon)
                                    latitudes.append(lat)
                                    datetimes.append(dt_str)
                                    categories.append(category_str)
                                except ValueError:
                                    continue
        
        if len(longitudes) == 0:
            print(f"Warning: No valid coordinates found in {kmz_file}")
            return np.array([]), np.array([]), [], []
            
        return np.array(longitudes), np.array(latitudes), datetimes, categories
        
    except zipfile.BadZipFile:
        print(f"Error: {kmz_file} is not a valid ZIP/KMZ file")
        return np.array([]), np.array([]), [], []
    except ET.ParseError as e:
        print(f"Error parsing KML file in {kmz_file}: {e}")
        return np.array([]), np.array([]), [], []
    except Exception as e:
        print(f"Error reading {kmz_file}: {e}")
        return np.array([]), np.array([]), [], []


def plot_solutions_2d(
        fig, ax,
        solution_file, timestep,
        variable='zeta', velocity_file=None,
        vmin=None, vmax=None, cmap='bwr',
        drawmesh=False, levels=20,
        title=None, plot_vectors=False,
        vector_scale=1.0, vector_color='black', vector_width=0.003,
        vector_legend=False, vector_legend_magnitude=None, vector_legend_location='southeast',
        vector_legend_unit='m/s', vector_legend_label='',
        xmin=None, xmax=None, ymin=None, ymax=None,
        cbar_increment=None, compute_disturbance=False, draw_shorelines=False,
        cbar_ticks_increment=None, cbar_label=None, vector_variable_x=None, vector_variable_y=None,
        vector_resample=False, vector_resample_dx=None, vector_resample_dy=None,
        track_file=None, track_color='red', track_linewidth=2.0, track_markersize=5.0,
        track_annotate_datetime=False, track_annotate_category=False, track_annotate_fontsize=8.0,
        compute_departure=False, departure_reference_file=None, departure_reference_time=None,
        departure_reference_variable=None,
        background_imagery=False, background_provider='esri', background_alpha=0.7,
        background_zoom=12, contour_alpha=1.0, show_axis_ticks=False, show_scale=False):
    """
    Plot 2D solution fields from CG ADCIRC NetCDF files.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to plot on
    ax : matplotlib.axes.Axes
        Axes object to plot on
    solution_file : str
        Path to the main solution NetCDF file (e.g., fort.63.nc)
    timestep : int
        Time step to plot (0-based index, use -1 for last time step)
    variable : str, optional
        Variable to plot (default: 'zeta'). Use 'velocity_mag' for velocity magnitude
    velocity_file : str, optional
        Path to velocity NetCDF file (fort.64.nc or fort.74.nc) for velocity plotting
    vector_variable_x : str, optional
        Name of variable for x-component of vectors in velocity_file (default: 'u-vel')
    vector_variable_y : str, optional
        Name of variable for y-component of vectors in velocity_file (default: 'v-vel')
    vmin : float, optional
        Minimum value for colorbar (auto if not specified)
    vmax : float, optional
        Maximum value for colorbar (auto if not specified)
    cmap : str, optional
        Colormap to use (default: 'bwr')
    drawmesh : bool, optional
        Draw mesh lines (default: False)
    levels : int, optional
        Number of contour levels (default: 20)
    title : str, optional
        Custom title for the plot
    plot_vectors : bool, optional
        Plot velocity vectors (default: False)
    vector_scale : float, optional
        Scale factor for velocity vectors (default: 1.0)
    vector_color : str, optional
        Color for velocity vectors (default: 'black')
    vector_width : float, optional
        Shaft width of velocity vector arrows (default: 0.003)
    vector_legend : bool, optional
        If True, draw a reference vector legend (default: False)
    vector_legend_magnitude : float, optional
        Magnitude of the reference vector in legend (required if vector_legend=True)
    vector_legend_location : str, optional
        Location for vector legend: 'northwest', 'northeast', 'southwest', 'southeast' (default: 'southeast')
    vector_legend_unit : str, optional
        Unit string to display in vector legend (default: 'm/s')
    vector_legend_label : str, optional
        Descriptive label text for the vector (e.g., 'Winds', 'Currents') (default: '')
    xmin : float, optional
        Minimum x value (longitude) for plot range
    xmax : float, optional
        Maximum x value (longitude) for plot range
    ymin : float, optional
        Minimum y value (latitude) for plot range
    ymax : float, optional
        Maximum y value (latitude) for plot range
    cbar_increment : float, optional
        Increment for colorbar levels (overrides levels parameter if specified)
    compute_disturbance : bool, optional
        If True, compute disturbance field on the fly as variable + min(depth, 0.0) (default: False)
    draw_shorelines : bool, optional
        If True, draw 0 m depth contour lines in gray (default: False)
    cbar_ticks_increment : float, optional
        Increment for colorbar tick labels (e.g., 0.5 creates ticks at vmin, vmin+0.5, ..., vmax)
    cbar_label : str, optional
        Custom label text for the colorbar (default: None, uses automatic labels based on variable)
    vector_variable_x : str, optional
        Name of variable for x-component of vectors in velocity_file (default: 'u-vel')
    vector_variable_y : str, optional
        Name of variable for y-component of vectors in velocity_file (default: 'v-vel')
    vector_resample : bool, optional
        If True, resample vectors at regular grid points (default: False)
    vector_resample_dx : float, optional
        Grid spacing in x direction for vector resampling (required if vector_resample=True)
    vector_resample_dy : float, optional
        Grid spacing in y direction for vector resampling (required if vector_resample=True)
    track_file : str, optional
        Path to KMZ file containing hurricane best track to overlay on plot
    track_color : str, optional
        Color for the track line (default: 'red')
    track_linewidth : float, optional
        Line width for the track (default: 2.0)
    track_markersize : float, optional
        Marker size for track points (default: 5.0)
    track_annotate_datetime : bool, optional
        Annotate track points with datetime labels (default: False)
    track_annotate_category : bool, optional
        Annotate track points with category labels (default: False)
    track_annotate_fontsize : float, optional
        Font size for track annotation text (default: 8.0)
    compute_departure : bool, optional
        If True, compute departure field by subtracting reference water level (default: False)
    departure_reference_file : str, optional
        Path to reference solution file (e.g., fort.63.nc) for departure computation
    departure_reference_time : str, optional
        Reference time in format "YYYY-MM-DD HH:mm:ss" (required if compute_departure=True)
    departure_reference_variable : str, optional
        Variable name to use from the reference file (default: 'zeta')
    background_imagery : bool, optional
        If True, add satellite/imagery background to the plot (default: False)
    background_provider : str, optional
        Background imagery provider: 'esri', 'osm', 'stamen', or 'google' (default: 'esri')
    background_alpha : float, optional
        Transparency of background imagery, 0.0 (transparent) to 1.0 (opaque) (default: 0.7)
    background_zoom : int, optional
        Zoom level for background imagery tiles. Higher values = higher resolution.
        Typical range: 8-15. Higher zoom levels may take longer to load (default: 12)
    contour_alpha : float, optional
        Transparency of the contour plot overlay, 0.0 (transparent) to 1.0 (opaque) (default: 1.0)
    show_axis_ticks : bool, optional
        If True, show axis ticks and labels when using background imagery (default: False)
    show_scale : bool, optional
        If True, show a scale bar on the plot (default: False)
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    
    # Check if files exist
    if not os.path.exists(solution_file):
        print(f"Error: File {solution_file} not found.")
        return False
    
    # Read the mesh data
    print(f"Reading mesh data from {solution_file}")
    try:
        ds = xr.open_dataset(solution_file)
    except Exception as e:
        print(f"Error reading {solution_file}: {e}")
        return False
    
    # Extract mesh information
    adc_x = ds.x.values
    adc_y = ds.y.values
    adc_e = ds.element.values - 1  # Convert to 0-based indexing
    adc_t = ds['time']  # Keep lazy - don't load all time values yet
    
    # Add background imagery if requested
    if background_imagery:
        if not CARTOPY_AVAILABLE:
            print("Warning: Cartopy not available. Background imagery disabled.")
        else:
            # Check if ax is a GeoAxes (should be if background_imagery is True)
            if not hasattr(ax, 'projection'):
                print("Warning: Background imagery requires GeoAxes. Creating one...")
                # This shouldn't happen if main() is called correctly, but handle it gracefully
                return False
            
            print(f"Adding background imagery from provider: {background_provider}")
            try:
                # Determine data extent for setting map bounds
                # Use user-specified limits if provided, otherwise use data bounds
                data_xmin = xmin if xmin is not None else np.min(adc_x)
                data_xmax = xmax if xmax is not None else np.max(adc_x)
                data_ymin = ymin if ymin is not None else np.min(adc_y)
                data_ymax = ymax if ymax is not None else np.max(adc_y)
                
                # Add padding to extent only if user didn't specify limits (10% on each side)
                if xmin is None and xmax is None:
                    x_pad = (data_xmax - data_xmin) * 0.1
                    data_xmin -= x_pad
                    data_xmax += x_pad
                if ymin is None and ymax is None:
                    y_pad = (data_ymax - data_ymin) * 0.1
                    data_ymin -= y_pad
                    data_ymax += y_pad
                
                extent = [data_xmin, data_xmax, data_ymin, data_ymax]
                
                # Add background imagery based on provider
                provider_lower = background_provider.lower()
                
                if provider_lower == 'esri':
                    # ESRI World Imagery using tile service (more reliable than WMS)
                    try:
                        # Create a custom tile source for ESRI World Imagery
                        # ESRI provides tile services at this URL pattern
                        class ESRIWorldImagery(cimgt.GoogleTiles):
                            def _image_url(self, tile):
                                x, y, z = tile
                                # ESRI World Imagery tile service
                                # Note: ESRI uses {z}/{y}/{x} format
                                url = ('https://server.arcgisonline.com/ArcGIS/rest/services/'
                                       'World_Imagery/MapServer/tile/{z}/{y}/{x}').format(
                                    z=z, y=y, x=x)
                                return url
                        
                        esri_tiles = ESRIWorldImagery()
                        ax.add_image(esri_tiles, background_zoom, alpha=background_alpha)
                    except Exception as e:
                        print(f"Warning: Could not load ESRI imagery: {e}")
                        print("Falling back to OSM imagery")
                        ax.add_image(OSM(), background_zoom, alpha=background_alpha)
                elif provider_lower == 'osm':
                    # OpenStreetMap
                    ax.add_image(OSM(), background_zoom, alpha=background_alpha)
                elif provider_lower == 'stamen':
                    # Stamen Terrain
                    from cartopy.io.img_tiles import StamenTerrain
                    ax.add_image(StamenTerrain(), background_zoom, alpha=background_alpha)
                elif provider_lower == 'google':
                    # Google (requires API key in some cases)
                    try:
                        from cartopy.io.img_tiles import GoogleTiles
                        ax.add_image(GoogleTiles(), background_zoom, alpha=background_alpha)
                    except Exception as e:
                        print(f"Warning: Could not load Google imagery: {e}")
                        print("Falling back to OSM imagery")
                        ax.add_image(OSM(), background_zoom, alpha=background_alpha)
                else:
                    print(f"Warning: Unknown background provider '{background_provider}'. Using OSM.")
                    ax.add_image(OSM(), background_zoom, alpha=background_alpha)
                
                # Set extent
                ax.set_extent(extent, crs=ccrs.PlateCarree())
                
            except Exception as e:
                print(f"Warning: Error adding background imagery: {e}")
                print("Continuing without background imagery...")
    
    # Handle -1 timestep (last time step)
    if timestep == -1:
        timestep_to_use = len(adc_t) - 1
        print(f"Using last time step: {timestep_to_use}")
    else:
        timestep_to_use = timestep
        # Check if the specified time step exists
        if timestep_to_use >= len(adc_t):
            print(f"Error: Time step {timestep_to_use} not found. Available time steps: 0 to {len(adc_t)-1}")
            ds.close()
            return False
    
    # Set default vector variable names if not specified
    if vector_variable_x is None:
        vector_variable_x = 'u-vel'
    if vector_variable_y is None:
        vector_variable_y = 'v-vel'
    
    # Handle special case for velocity_mag
    if variable == 'velocity_mag':
        if velocity_file is None:
            print("Error: velocity_file must be specified for velocity_mag plotting")
            ds.close()
            return False
            
        if not os.path.exists(velocity_file):
            print(f"Error: Velocity file {velocity_file} not found.")
            ds.close()
            return False
        
        print(f"Reading velocity data from {velocity_file}")
        try:
            ds_vel = xr.open_dataset(velocity_file)
        except Exception as e:
            print(f"Error reading {velocity_file}: {e}")
            ds.close()
            return False
        
        # Check if velocity variables exist
        if vector_variable_x not in ds_vel.variables or vector_variable_y not in ds_vel.variables:
            print(f"Error: {vector_variable_x} or {vector_variable_y} not found in {velocity_file}")
            print(f"Available variables: {list(ds_vel.variables.keys())}")
            ds_vel.close()
            ds.close()
            return False
        
        # Extract velocity components (index first, then call .values for efficiency)
        u_vel = ds_vel[vector_variable_x][timestep_to_use, :].values
        v_vel = ds_vel[vector_variable_y][timestep_to_use, :].values
        
        # Calculate velocity magnitude
        var_data = np.sqrt(u_vel**2 + v_vel**2)
        
        # Close velocity dataset
        ds_vel.close()
        
    else:
        # Check if the specified variable exists
        if variable not in ds.variables:
            print(f"Error: Variable '{variable}' not found in {solution_file}")
            print(f"Available variables: {list(ds.variables.keys())}")
            ds.close()
            return False
        
        # Extract the variable data for the specified time step (index first, then call .values for efficiency)
        var_shape = ds[variable].shape
        if len(var_shape) == 2:
            # Variables with shape (time, node) - defined on nodes
            var_data = ds[variable][timestep_to_use, :].values
        elif len(var_shape) == 3 and var_shape[2] == 1:
            # Variables with shape (time, element, 1) - defined on elements
            var_data = ds[variable][timestep_to_use, :, 0].values
        else:
            # For other cases, try to get the data at the specified time step
            var_data = ds[variable][timestep_to_use, :].values
    
    # Handle velocity data for vector plotting
    if plot_vectors:
        if velocity_file is None:
            print("Error: velocity_file must be specified for vector plotting")
            ds.close()
            return False
            
        if not os.path.exists(velocity_file):
            print(f"Error: Velocity file {velocity_file} not found for vector plotting.")
            ds.close()
            return False
        
        print(f"Reading velocity data from {velocity_file}")
        try:
            ds_vel = xr.open_dataset(velocity_file)
        except Exception as e:
            print(f"Error reading {velocity_file}: {e}")
            ds.close()
            return False
        
        # Check if velocity variables exist
        if vector_variable_x not in ds_vel.variables or vector_variable_y not in ds_vel.variables:
            print(f"Error: {vector_variable_x} or {vector_variable_y} not found in {velocity_file}")
            print(f"Available variables: {list(ds_vel.variables.keys())}")
            ds_vel.close()
            ds.close()
            return False
        
        # Extract velocity components (index first, then call .values for efficiency)
        u_vel = ds_vel[vector_variable_x][timestep_to_use, :].values
        v_vel = ds_vel[vector_variable_y][timestep_to_use, :].values
        
        # Close velocity dataset
        ds_vel.close()
    
    # Read depth data if needed for disturbance or shorelines
    depth = None
    if compute_disturbance or draw_shorelines:
        # Check if depth variable exists
        if 'depth' not in ds.variables:
            if compute_disturbance:
                print(f"Error: 'depth' variable not found in {solution_file} (required for disturbance computation)")
                ds.close()
                return False
            elif draw_shorelines:
                print(f"Error: 'depth' variable not found in {solution_file} (required for shoreline drawing)")
                ds.close()
                return False
        
        # Extract depth data (should be on nodes, no time dimension)
        depth = ds['depth'].values
    
    # Compute disturbance if requested
    if compute_disturbance:
        # Disturbance computation only works for variables from the solution file, not velocity_mag
        if variable == 'velocity_mag':
            print(f"Error: Cannot compute disturbance for 'velocity_mag'. Use a variable from the solution file (e.g., 'zeta').")
            ds.close()
            return False
        
        # Check dimensions match
        if len(depth) != len(var_data):
            print(f"Error: Dimension mismatch between '{variable}' and 'depth' for disturbance computation")
            ds.close()
            return False
        
        # Compute disturbance: variable + min(depth, 0.0)
        print(f"Computing disturbance field from '{variable}' and 'depth'")
        disturbance = var_data + np.minimum(depth, 0.0)
        var_data = disturbance
        # Update variable name for plotting
        variable = 'disturbance'
    
    # Compute departure if requested
    if compute_departure:
        # Departure computation only works for variables from the solution file, not velocity_mag
        if variable == 'velocity_mag':
            print(f"Error: Cannot compute departure for 'velocity_mag'. Use a variable from the solution file (e.g., 'zeta').")
            ds.close()
            return False
        
        if departure_reference_file is None or departure_reference_time is None:
            print(f"Error: departure_reference_file and departure_reference_time must be specified when compute_departure=True")
            ds.close()
            return False
        
        if not os.path.exists(departure_reference_file):
            print(f"Error: Reference file {departure_reference_file} not found.")
            ds.close()
            return False
        
        # Check if depth is available (needed for dry node handling)
        if depth is None:
            if 'depth' not in ds.variables:
                print(f"Error: 'depth' variable not found in {solution_file} (required for departure computation)")
                ds.close()
                return False
            depth = ds['depth'].values
        
        print(f"Reading reference data from {departure_reference_file}")
        try:
            ds_ref = xr.open_dataset(departure_reference_file)
        except Exception as e:
            print(f"Error reading {departure_reference_file}: {e}")
            ds.close()
            return False
        
        # Determine which variable to use from reference file
        ref_variable = departure_reference_variable if departure_reference_variable is not None else 'zeta'
        
        # Check if reference variable exists
        if ref_variable not in ds_ref.variables:
            print(f"Error: Variable '{ref_variable}' not found in {departure_reference_file}")
            print(f"Available variables: {list(ds_ref.variables.keys())}")
            ds_ref.close()
            ds.close()
            return False
        
        # Parse reference time
        try:
            # Try parsing with different formats
            try:
                ref_time = pd.to_datetime(departure_reference_time)
            except:
                # Try datetime parsing
                ref_time = datetime.strptime(departure_reference_time, "%Y-%m-%d %H:%M:%S")
        except Exception as e:
            print(f"Error: Cannot parse departure_reference_time '{departure_reference_time}'. Use format 'YYYY-MM-DD HH:mm:ss'")
            ds_ref.close()
            ds.close()
            return False
        
        # Get time values from reference dataset
        ref_time_values = ds_ref['time'].values
        
        # Convert to pandas Timestamp for comparison
        if isinstance(ref_time_values[0], np.datetime64):
            ref_time_pd = pd.to_datetime(ref_time)
            # Find closest time step using vectorized operations
            ref_time_values_pd = pd.to_datetime(ref_time_values)
            time_diffs = np.abs(ref_time_values_pd - ref_time_pd)
            ref_timestep = np.argmin(time_diffs)
        else:
            # Try to convert and find matching time
            print(f"Warning: Unexpected time format in reference file. Attempting to find matching time step.")
            ref_timestep = 0
        
        print(f"Using reference time step {ref_timestep} (target: {departure_reference_time})")
        print(f"Using reference variable '{ref_variable}' from reference file")
        
        # Extract reference variable data (index first, then call .values for efficiency)
        ref_shape = ds_ref[ref_variable].shape
        if len(ref_shape) == 2:
            zeta_ref = ds_ref[ref_variable][ref_timestep, :].values
        else:
            print(f"Error: Unexpected shape for '{ref_variable}' in reference file: {ref_shape}")
            ds_ref.close()
            ds.close()
            return False
        
        # Check dimensions match
        if len(zeta_ref) != len(var_data) or len(zeta_ref) != len(depth):
            print(f"Error: Dimension mismatch for departure computation")
            print(f"  Current variable: {len(var_data)}, Reference: {len(zeta_ref)}, Depth: {len(depth)}")
            ds_ref.close()
            ds.close()
            return False
        
        # Close reference dataset
        ds_ref.close()
        
        # Compute reference values: use zeta_ref for wet nodes, use elevation (-depth) for dry nodes
        # A node is dry if depth < 0
        reference_values = np.where(depth < 0, -depth, zeta_ref)
        
        # Compute departure: current - reference
        print(f"Computing departure field from '{variable}' and reference values")
        departure = var_data - reference_values
        var_data = departure
        # Update variable name for plotting
        variable = 'departure'
    
    # Close the main dataset
    ds.close()
    
    # Set up colorbar limits
    vmin_plot = vmin if vmin is not None else np.nanmin(var_data)
    vmax_plot = vmax if vmax is not None else np.nanmax(var_data)
    
    # Create contour plot
    # Check for non-finite values
    finite_mask = np.isfinite(var_data)
    if not np.any(finite_mask):
        print("Warning: All data values are non-finite. Cannot create plot.")
        return False
    
    # Create triangulation and mask triangles with non-finite values
    triang = mtri.Triangulation(adc_x, adc_y, adc_e)
    
    # Mask triangles that contain non-finite values
    if len(var_data) == len(adc_x):  # Nodal data
        # Mask triangles where any vertex has non-finite value (reuse finite_mask)
        mask = np.any(~finite_mask[adc_e], axis=1)
    else:  # Elemental data
        # For elemental data, mask elements with non-finite values (reuse finite_mask)
        mask = ~finite_mask
    
    triang.set_mask(mask)
    
    # For nodal data, use tricontourf
    if len(var_data) == len(adc_x):  # Nodal data
        # Create contour plot with levels
        if cbar_increment is not None:
            # Use specified increment for levels
            levels_plot = np.arange(vmin_plot, vmax_plot + cbar_increment, cbar_increment)
        else:
            # Use specified number of levels
            levels_plot = np.linspace(vmin_plot, vmax_plot, levels)
        contour = ax.tricontourf(triang, var_data, levels=levels_plot, cmap=cmap, extend='both', alpha=contour_alpha)
    else:  # Elemental data
        contour = ax.tripcolor(triang, var_data, cmap=cmap, vmin=vmin_plot, vmax=vmax_plot, alpha=contour_alpha)
    
    # Draw mesh if requested
    if drawmesh:
        ax.triplot(triang, color='black', linewidth=0.1, alpha=0.5)
    
    # Draw shorelines (0 m depth contour) if requested
    if draw_shorelines:
        if depth is None:
            print("Warning: Depth data not available for shoreline drawing. Skipping shorelines.")
        else:
            # Create a separate triangulation for depth (not masked by var_data)
            # Mask triangles where depth is non-finite for shoreline drawing
            if len(depth) == len(adc_x):  # Nodal depth data
                depth_mask = np.any(~np.isfinite(depth[adc_e]), axis=1)
            else:  # Elemental depth data
                depth_mask = ~np.isfinite(depth)
            
            depth_triang = mtri.Triangulation(adc_x, adc_y, adc_e)
            depth_triang.set_mask(depth_mask)
            
            # Draw contour line at depth = 0
            ax.tricontour(depth_triang, depth, levels=[0.0], colors='black', linewidths=0.1, linestyles='-', zorder=5)
    
    # Plot velocity vectors if requested
    quiver_plot = None  # Store quiver plot for legend
    if plot_vectors:
        # Check if velocity data is available
        if 'u_vel' not in locals() or 'v_vel' not in locals():
            print("Error: Velocity data not available for vector plotting.")
            return False
        
        # Resample vectors if requested
        if vector_resample:
            if vector_resample_dx is None or vector_resample_dy is None:
                print("Error: vector_resample_dx and vector_resample_dy must be specified when vector_resample=True")
                return False
            
            # Determine grid extent (use plot limits if specified, otherwise use data extent)
            x_min = xmin if xmin is not None else np.min(adc_x)
            x_max = xmax if xmax is not None else np.max(adc_x)
            y_min = ymin if ymin is not None else np.min(adc_y)
            y_max = ymax if ymax is not None else np.max(adc_y)
            
            # Create regular grid, shifted by half dx and dy to avoid vectors at boundaries
            x_grid = np.arange(x_min + vector_resample_dx/2, x_max, vector_resample_dx)
            y_grid = np.arange(y_min + vector_resample_dy/2, y_max, vector_resample_dy)
            X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
            
            # Prepare input points (only use finite values for interpolation)
            valid_mask = np.isfinite(u_vel) & np.isfinite(v_vel)
            points = np.column_stack((adc_x[valid_mask], adc_y[valid_mask]))
            u_values = u_vel[valid_mask]
            v_values = v_vel[valid_mask]
            
            # Interpolate u and v components to regular grid
            u_grid = griddata(points, u_values, (X_grid, Y_grid), method='linear', fill_value=np.nan)
            v_grid = griddata(points, v_values, (X_grid, Y_grid), method='linear', fill_value=np.nan)
            
            # Mask out NaN values
            valid_grid = np.isfinite(u_grid) & np.isfinite(v_grid)
            
            # Plot resampled vectors
            quiver_plot = ax.quiver(X_grid[valid_grid], Y_grid[valid_grid], 
                     u_grid[valid_grid], v_grid[valid_grid],
                     scale=1.0/vector_scale, scale_units='xy', 
                     color=vector_color, width=vector_width, alpha=0.8)
        else:
            # Plot vectors at original nodes
            quiver_plot = ax.quiver(adc_x, adc_y, u_vel, v_vel, 
                     scale=1.0/vector_scale, scale_units='xy', 
                     color=vector_color, width=vector_width, alpha=0.8)
        
        # Add vector legend if requested
        if vector_legend:
            if vector_legend_magnitude is None:
                print("Warning: vector_legend_magnitude must be specified when vector_legend=True. Skipping legend.")
            else:
                # Map location string to coordinates
                location_map = {
                    'northwest': (0.1, 0.95),
                    'northeast': (0.9, 0.95),
                    'southwest': (0.1, 0.05),
                    'southeast': (0.9, 0.05)
                }
                
                if vector_legend_location not in location_map:
                    print(f"Warning: Invalid vector_legend_location '{vector_legend_location}'. Using 'southeast'.")
                    vector_legend_location = 'southeast'
                
                loc_x, loc_y = location_map[vector_legend_location]
                
                # Construct legend label
                if vector_legend_label:
                    legend_text = f'{vector_legend_label}: {vector_legend_magnitude} {vector_legend_unit}'
                else:
                    legend_text = f'{vector_legend_magnitude} {vector_legend_unit}'
                
                # Add quiver key (reference arrow)
                # labelpos='E' places the label to the right of the arrow
                qk = ax.quiverkey(quiver_plot, loc_x, loc_y, vector_legend_magnitude,
                            legend_text,
                            labelpos='E', coordinates='axes',
                            color=vector_color, fontproperties={'size': 8})
    
    # Plot hurricane track if requested
    if track_file:
        if os.path.exists(track_file):
            print(f"Reading hurricane track from {track_file}")
            track_lons, track_lats, track_dts, track_categories = parse_kmz_track(track_file)
            if len(track_lons) > 0:
                # Plot track line
                ax.plot(track_lons, track_lats, color=track_color, 
                       linewidth=track_linewidth, linestyle='-', 
                       marker='o', markersize=track_markersize, 
                       label='Hurricane Track', zorder=10)
                
                # Add annotations if requested
                if track_annotate_datetime or track_annotate_category:
                    for lon, lat, dt, cat in zip(track_lons, track_lats, track_dts, track_categories):
                        annotation_parts = []
                        
                        # Add datetime if requested
                        if track_annotate_datetime and dt:
                            # Format datetime to short format
                            dt_formatted = format_datetime_short(dt)
                            annotation_parts.append(dt_formatted)
                        
                        # Add category if requested
                        if track_annotate_category:
                            category_suffix = get_category_suffix(cat)
                            if category_suffix:
                                # Remove leading comma and space if datetime was not added
                                if not track_annotate_datetime:
                                    category_suffix = category_suffix.lstrip(', ')
                                annotation_parts.append(category_suffix)
                        
                        # Create annotation text
                        if annotation_parts:
                            # Join parts: use empty string to preserve spacing in category_suffix (which includes comma)
                            annotation_text = ''.join(annotation_parts)
                            ax.annotate(annotation_text, (lon, lat), xytext=(5, 5), 
                                      textcoords='offset points', fontsize=track_annotate_fontsize,
                                      zorder=11)
            else:
                print(f"Warning: No track data extracted from {track_file}")
        else:
            print(f"Warning: Track file {track_file} not found")
    
    # Add colorbar
    cbar = fig.colorbar(contour, ax=ax)
    if cbar_label:
        cbar.ax.set_ylabel(cbar_label, rotation=270, labelpad=15)
    elif variable == 'velocity_mag':
        cbar.ax.set_ylabel('Velocity Magnitude (m/s)', rotation=270, labelpad=15)
    elif variable == 'zeta':
        cbar.ax.set_ylabel('Water Level (m)', rotation=270, labelpad=15)
    elif variable == 'disturbance':
        cbar.ax.set_ylabel('Disturbance (m)', rotation=270, labelpad=15)
    elif variable == 'departure':
        cbar.ax.set_ylabel('Departure (m)', rotation=270, labelpad=15)
    else:
        cbar.ax.set_ylabel(f'{variable} (units)', rotation=270, labelpad=15)
    
    # Set colorbar tick increment if specified
    if cbar_ticks_increment is not None:
        # Generate ticks from vmin_plot to vmax_plot with specified increment
        ticks = np.arange(vmin_plot, vmax_plot + cbar_ticks_increment, cbar_ticks_increment)
        # Ensure vmax_plot is included if it's close to the last tick
        if len(ticks) > 0 and abs(ticks[-1] - vmax_plot) > cbar_ticks_increment * 0.01:
            ticks = np.append(ticks, vmax_plot)
        # Always include 0.0 if it's within [vmin_plot, vmax_plot]
        if vmin_plot <= 0.0 <= vmax_plot:
            # Check if 0.0 is already in ticks (within tolerance)
            if len(ticks) == 0 or min(np.abs(ticks - 0.0)) > cbar_ticks_increment * 0.01:
                ticks = np.append(ticks, 0.0)
                ticks = np.sort(ticks)
        cbar.set_ticks(ticks)
    
    # Set labels and title
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    
    if title:
        ax.set_title(title)
    else:
        time_str = str(adc_t[timestep_to_use].values)[:19].replace('T', ' ') + ' UTC'
        ax.set_title(time_str)
    
    # Set plot limits if specified
    # For GeoAxes, use set_extent; for regular axes, use set_xlim/set_ylim
    is_geoaxes = hasattr(ax, 'projection')
    
    if is_geoaxes:
        # GeoAxes: use set_extent with PlateCarree CRS
        if xmin is not None or xmax is not None or ymin is not None or ymax is not None:
            # Get current extent if limits not specified
            current_extent = ax.get_extent(crs=ccrs.PlateCarree())
            xlim = [xmin if xmin is not None else current_extent[0],
                    xmax if xmax is not None else current_extent[1]]
            ylim = [ymin if ymin is not None else current_extent[2],
                    ymax if ymax is not None else current_extent[3]]
            ax.set_extent([xlim[0], xlim[1], ylim[0], ylim[1]], crs=ccrs.PlateCarree())
        # Aspect ratio is handled automatically by GeoAxes based on projection
        
        # Enable gridlines with tick labels for GeoAxes if requested
        if show_axis_ticks:
            gl = ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,
                             linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': 10}
            gl.ylabel_style = {'size': 10}
        
        # Add scale bar if requested (for GeoAxes)
        if show_scale:
            try:
                if SCALEBAR_AVAILABLE:
                    # Use matplotlib-scalebar library
                    # Calculate approximate scale based on extent
                    extent = ax.get_extent(crs=ccrs.PlateCarree())
                    center_lat = (extent[2] + extent[3]) / 2
                    meters_per_degree = 111320 * np.cos(np.radians(center_lat))
                    
                    # Choose a reasonable length (10% of map width)
                    map_width_deg = extent[1] - extent[0]
                    scale_length_m = map_width_deg * meters_per_degree * 0.1
                    
                    # Round to nice value
                    if scale_length_m < 1000:
                        scale_length = round(scale_length_m / 100) * 100  # Round to 100m
                        unit = 'm'
                    elif scale_length_m < 10000:
                        scale_length = round(scale_length_m / 1000)  # Round to km
                        unit = 'km'
                    else:
                        scale_length = round(scale_length_m / 10000) * 10  # Round to 10km
                        unit = 'km'
                    
                    # Convert to meters for ScaleBar
                    scale_length_meters = scale_length * (1000 if unit == 'km' else 1)
                    
                    # Use dx parameter for ScaleBar (in meters per pixel, but we'll use approximate)
                    # For geographic coordinates, we need to specify dx in the data units
                    # Since ScaleBar works with data coordinates, we'll use a workaround
                    # by specifying the length in the projection's units
                    scalebar = ScaleBar(scale_length_meters, units='m', location='lower left',
                                      box_alpha=0.7, font_properties={'size': 10},
                                      length_fraction=0.2, width_fraction=0.01,
                                      dx=1.0, frameon=True)
                    # Note: ScaleBar with geographic coordinates is tricky, so we'll use manual fallback
                    # Actually, let's use the manual implementation which is more reliable for GeoAxes
                    raise ImportError("Using manual scale bar for GeoAxes")
                else:
                    raise ImportError("ScaleBar library not available")
            except (ImportError, AttributeError):
                # Manual scale bar implementation for GeoAxes
                extent = ax.get_extent(crs=ccrs.PlateCarree())
                center_lat = (extent[2] + extent[3]) / 2
                meters_per_degree = 111320 * np.cos(np.radians(center_lat))
                
                # Determine scale bar length and position
                map_width_deg = extent[1] - extent[0]
                scale_length_deg = map_width_deg * 0.1  # 10% of width
                scale_length_m = scale_length_deg * meters_per_degree
                
                # Round to nice value
                if scale_length_m < 1000:
                    scale_length_m = round(scale_length_m / 100) * 100
                    label = f"{int(scale_length_m)} m"
                else:
                    scale_length_km = round(scale_length_m / 1000)
                    scale_length_m = scale_length_km * 1000
                    label = f"{scale_length_km} km"
                
                scale_length_deg = scale_length_m / meters_per_degree
                
                # Position at lower left (with some padding)
                x_pos = extent[0] + (extent[1] - extent[0]) * 0.05
                y_pos = extent[2] + (extent[3] - extent[2]) * 0.05
                
                # Draw scale bar
                ax.plot([x_pos, x_pos + scale_length_deg], [y_pos, y_pos], 
                       'k-', linewidth=3, transform=ccrs.PlateCarree(), zorder=20)
                ax.text(x_pos + scale_length_deg/2, y_pos + (extent[3] - extent[2]) * 0.02,
                       label, ha='center', va='bottom', transform=ccrs.PlateCarree(),
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7),
                       fontsize=10, zorder=20)
            except Exception as e:
                print(f"Warning: Could not add scale bar: {e}")
    else:
        # Regular axes: use set_xlim/set_ylim
        if xmin is not None or xmax is not None:
            xlim = [xmin if xmin is not None else ax.get_xlim()[0],
                    xmax if xmax is not None else ax.get_xlim()[1]]
            ax.set_xlim(xlim)
        
        if ymin is not None or ymax is not None:
            ylim = [ymin if ymin is not None else ax.get_ylim()[0],
                    ymax if ymax is not None else ax.get_ylim()[1]]
            ax.set_ylim(ylim)
        
        # Set aspect ratio for regular axes
        ax.set_aspect('equal')
        
        # Add scale bar if requested (for regular axes)
        if show_scale:
            try:
                if SCALEBAR_AVAILABLE:
                    # For regular axes, assume coordinates are in meters
                    # Calculate scale based on data extent
                    x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
                    scale_length = x_range * 0.1  # 10% of plot width
                    
                    if scale_length < 1000:
                        scale_length = round(scale_length / 100) * 100
                        unit = 'm'
                    else:
                        scale_length = round(scale_length / 1000)
                        unit = 'km'
                    
                    scalebar = ScaleBar(scale_length, unit=unit, location='lower left',
                                      box_alpha=0.7, font_properties={'size': 10},
                                      length_fraction=0.2, width_fraction=0.01)
                    ax.add_artist(scalebar)
                else:
                    # Manual scale bar for regular axes
                    x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
                    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
                    scale_length = x_range * 0.1
                    
                    if scale_length < 1000:
                        scale_length = round(scale_length / 100) * 100
                        label = f"{int(scale_length)} m"
                    else:
                        scale_length_km = round(scale_length / 1000)
                        scale_length = scale_length_km * 1000
                        label = f"{scale_length_km} km"
                    
                    # Position at lower left
                    x_pos = ax.get_xlim()[0] + x_range * 0.05
                    y_pos = ax.get_ylim()[0] + y_range * 0.05
                    
                    # Draw scale bar
                    ax.plot([x_pos, x_pos + scale_length], [y_pos, y_pos], 
                           'k-', linewidth=3, zorder=20)
                    ax.text(x_pos + scale_length/2, y_pos + y_range * 0.02,
                           label, ha='center', va='bottom',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7),
                           fontsize=10, zorder=20)
            except Exception as e:
                print(f"Warning: Could not add scale bar: {e}")
    
    return True


def get_parser(add_help=True):
    """Get command line argument parser."""
    parser = argparse.ArgumentParser(description='Plot CG ADCIRC water levels and velocity fields from NetCDF files.', add_help=add_help)
    parser.add_argument('--solution-file', type=str, required=True, help='Path to main solution NetCDF file (e.g., fort.63.nc)')
    parser.add_argument('--timestep', type=int, required=True, help='Time step to plot (0-based index, use -1 for last time step)')
    parser.add_argument('--variable', type=str, default='zeta', help='Variable to plot (default: zeta). Use "velocity_mag" for velocity magnitude')
    parser.add_argument('--velocity-file', type=str, help='Path to velocity NetCDF file (fort.64.nc or fort.74.nc) for velocity plotting')
    parser.add_argument('--output', type=str, required=True, help='Output PNG file name')
    parser.add_argument('--figsizex', type=float, default=12.0, help='Figure size in x direction')
    parser.add_argument('--figsizey', type=float, default=10.0, help='Figure size in y direction')
    parser.add_argument('--vmin', type=float, help='Minimum value of the colorbar (auto if not specified)')
    parser.add_argument('--vmax', type=float, help='Maximum value of the colorbar (auto if not specified)')
    parser.add_argument('--cmap', type=str, default='bwr', help='Colormap to use (default: bwr)')
    parser.add_argument('--drawmesh', action='store_true', help='Draw mesh lines')
    parser.add_argument('--levels', type=int, default=20, help='Number of contour levels (default: 20)')
    parser.add_argument('--title', type=str, help='Custom title for the plot')
    parser.add_argument('--xmin', type=float, help='Minimum x value (longitude) for plot range')
    parser.add_argument('--xmax', type=float, help='Maximum x value (longitude) for plot range')
    parser.add_argument('--ymin', type=float, help='Minimum y value (latitude) for plot range')
    parser.add_argument('--ymax', type=float, help='Maximum y value (latitude) for plot range')
    parser.add_argument('--plot-vectors', action='store_true', help='Plot velocity vectors')
    parser.add_argument('--vector-scale', type=float, default=1.0, help='Scale factor for velocity vectors (default: 1.0)')
    parser.add_argument('--vector-color', type=str, default='black', help='Color for velocity vectors (default: black)')
    parser.add_argument('--vector-width', type=float, default=0.003, help='Shaft width of velocity vector arrows (default: 0.003)')
    parser.add_argument('--vector-legend', action='store_true', help='Draw a reference vector legend')
    parser.add_argument('--vector-legend-magnitude', type=float, help='Magnitude of the reference vector in legend')
    parser.add_argument('--vector-legend-location', type=str, default='southeast', 
                       choices=['northwest', 'northeast', 'southwest', 'southeast'],
                       help='Location for vector legend (default: southeast)')
    parser.add_argument('--vector-legend-unit', type=str, default='m/s', help='Unit string to display in vector legend (default: m/s)')
    parser.add_argument('--vector-legend-label', type=str, default='', help='Descriptive label text for the vector (e.g., Winds, Currents) (default: blank)')
    parser.add_argument('--vector-variable-x', type=str, help='Name of variable for x-component of vectors in velocity_file (default: u-vel)')
    parser.add_argument('--vector-variable-y', type=str, help='Name of variable for y-component of vectors in velocity_file (default: v-vel)')
    parser.add_argument('--vector-resample', action='store_true', help='Resample vectors at regular grid points for plotting')
    parser.add_argument('--vector-resample-dx', type=float, help='Grid spacing in x direction for vector resampling (required if --vector-resample is used)')
    parser.add_argument('--vector-resample-dy', type=float, help='Grid spacing in y direction for vector resampling (required if --vector-resample is used)')
    parser.add_argument('--cbar-increment', type=float, help='Increment for colorbar levels (overrides --levels if specified)')
    parser.add_argument('--cbar-ticks-increment', type=float, help='Increment for colorbar tick labels (e.g., 0.5 creates ticks at vmin, vmin+0.5, ..., vmax)')
    parser.add_argument('--cbar-label', type=str, help='Custom label text for the colorbar (default: automatic based on variable)')
    parser.add_argument('--disturbance', action='store_true', help='Compute and plot disturbance field on the fly (disturbance = variable + min(depth, 0.0))')
    parser.add_argument('--departure', action='store_true', help='Compute and plot departure field by subtracting reference water level')
    parser.add_argument('--departure-reference-file', type=str, help='Path to reference solution file (e.g., fort.63.nc) for departure computation')
    parser.add_argument('--departure-reference-time', type=str, help='Reference time in format "YYYY-MM-DD HH:mm:ss" (required if --departure is used)')
    parser.add_argument('--departure-reference-variable', type=str, default='zeta', help='Variable name to use from the reference file (default: zeta)')
    parser.add_argument('--draw-shorelines', action='store_true', help='Draw 0 m depth contour lines in gray')
    parser.add_argument('--track-file', type=str, help='Path to KMZ file containing hurricane best track to overlay on plot')
    parser.add_argument('--track-color', type=str, default='red', help='Color for the track line (default: red)')
    parser.add_argument('--track-linewidth', type=float, default=2.0, help='Line width for the track (default: 2.0)')
    parser.add_argument('--track-markersize', type=float, default=5.0, help='Marker size for track points (default: 5.0)')
    parser.add_argument('--track-annotate-datetime', action='store_true', help='Annotate track points with datetime labels')
    parser.add_argument('--track-annotate-category', action='store_true', help='Annotate track points with category labels')
    parser.add_argument('--track-annotate-fontsize', type=float, default=8.0, help='Font size for track annotation text (default: 8.0)')
    parser.add_argument('--background-imagery', action='store_true', help='Add satellite/imagery background to the plot')
    parser.add_argument('--background-provider', type=str, default='esri', 
                       choices=['esri', 'osm', 'stamen', 'google'],
                       help='Background imagery provider: esri, osm, stamen, or google (default: esri)')
    parser.add_argument('--background-alpha', type=float, default=0.7, 
                       help='Transparency of background imagery, 0.0 (transparent) to 1.0 (opaque) (default: 0.7)')
    parser.add_argument('--background-zoom', type=int, default=12,
                       help='Zoom level for background imagery tiles. Higher values = higher resolution. Typical range: 8-15. Higher zoom levels may take longer to load (default: 12)')
    parser.add_argument('--contour-alpha', type=float, default=1.0,
                       help='Transparency of the contour plot overlay, 0.0 (transparent) to 1.0 (opaque) (default: 1.0)')
    parser.add_argument('--show-axis-ticks', action='store_true',
                       help='Show axis ticks and labels when using background imagery (default: False)')
    parser.add_argument('--show-scale', action='store_true',
                       help='Show a scale bar on the plot (default: False)')
    
    return parser


def main(args=None):
    """Main function for command line usage."""
    if args is None:
        args = get_parser().parse_args()
    
    # Validate velocity file requirement for certain options
    if args.variable == 'velocity_mag' and args.velocity_file is None:
        print("Error: --velocity-file must be specified when using --variable velocity_mag")
        sys.exit(1)
    
    if args.plot_vectors and args.velocity_file is None:
        print("Error: --velocity-file must be specified when using --plot-vectors")
        sys.exit(1)
    
    # Validate vector resampling options
    if args.vector_resample:
        if args.vector_resample_dx is None or args.vector_resample_dy is None:
            print("Error: --vector-resample-dx and --vector-resample-dy must be specified when using --vector-resample")
            sys.exit(1)
        if not args.plot_vectors:
            print("Warning: --vector-resample is ignored when --plot-vectors is not used")
    
    # Validate vector legend options
    if args.vector_legend:
        if not args.plot_vectors:
            print("Error: --vector-legend requires --plot-vectors")
            sys.exit(1)
        if args.vector_legend_magnitude is None:
            print("Error: --vector-legend-magnitude must be specified when using --vector-legend")
            sys.exit(1)
    
    # Validate departure options
    if args.departure:
        if args.departure_reference_file is None or args.departure_reference_time is None:
            print("Error: --departure-reference-file and --departure-reference-time must be specified when using --departure")
            sys.exit(1)
        if args.variable == 'velocity_mag':
            print("Error: Cannot compute departure for 'velocity_mag'. Use a variable from the solution file (e.g., 'zeta').")
            sys.exit(1)
    
    # Create the plot
    print(f"Creating plot for time step {args.timestep}")
    
    # Create GeoAxes if background imagery is enabled, otherwise regular axes
    if args.background_imagery:
        if not CARTOPY_AVAILABLE:
            print("Error: Cartopy is required for background imagery but is not available.")
            print("Please install cartopy: pip install cartopy")
            sys.exit(1)
        fig = plt.figure(figsize=(args.figsizex, args.figsizey))
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    else:
        fig, ax = plt.subplots(figsize=(args.figsizex, args.figsizey))
    
    success = plot_solutions_2d(
        fig, ax,
        args.solution_file, args.timestep,
        variable=args.variable, velocity_file=args.velocity_file,
        vmin=args.vmin, vmax=args.vmax, cmap=args.cmap,
        drawmesh=args.drawmesh, levels=args.levels,
        title=args.title, plot_vectors=args.plot_vectors,
        vector_scale=args.vector_scale, vector_color=args.vector_color,
        vector_width=args.vector_width,
        vector_legend=args.vector_legend, vector_legend_magnitude=args.vector_legend_magnitude,
        vector_legend_location=args.vector_legend_location, vector_legend_unit=args.vector_legend_unit,
        vector_legend_label=args.vector_legend_label,
        xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax,
        cbar_increment=args.cbar_increment, compute_disturbance=args.disturbance,
        draw_shorelines=args.draw_shorelines, cbar_ticks_increment=args.cbar_ticks_increment,
        cbar_label=args.cbar_label, vector_variable_x=args.vector_variable_x, vector_variable_y=args.vector_variable_y,
        vector_resample=args.vector_resample, vector_resample_dx=args.vector_resample_dx,
        vector_resample_dy=args.vector_resample_dy,
        track_file=args.track_file, track_color=args.track_color,
        track_linewidth=args.track_linewidth, track_markersize=args.track_markersize,
        track_annotate_datetime=args.track_annotate_datetime,
        track_annotate_category=args.track_annotate_category,
        track_annotate_fontsize=args.track_annotate_fontsize,
        compute_departure=args.departure, departure_reference_file=args.departure_reference_file,
        departure_reference_time=args.departure_reference_time,
        departure_reference_variable=args.departure_reference_variable,
        background_imagery=args.background_imagery, background_provider=args.background_provider,
        background_alpha=args.background_alpha, background_zoom=args.background_zoom,
        contour_alpha=args.contour_alpha, show_axis_ticks=args.show_axis_ticks,
        show_scale=args.show_scale
    )
    
    if not success:
        print("Error: Failed to create plot")
        sys.exit(1)
    
    # Save the plot
    print(f"Saving plot to {args.output}")
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print("Done!")


if __name__ == "__main__":
    main()
