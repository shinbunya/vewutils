#!/usr/bin/env python3
"""
Plot CG ADCIRC maximum water level fields from maxele.63.nc NetCDF files.
Similar functionality to plot_solution_2d.py but for maximum water level data.
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
from dateutil import parser as date_parser
import zipfile
import tempfile
import xml.etree.ElementTree as ET


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
        from dateutil import parser as date_parser
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


def get_category_number(category_str, include_non_hurricane=False):
    """
    Convert category string to compact label for inside-circle annotation.

    Parameters
    ----------
    category_str : str
        Category string (e.g., 'cat1', 'cat2', 'ts', 'ex')
    include_non_hurricane : bool, optional
        If True, include non-hurricane categories such as TS and EX (default: False)

    Returns
    -------
    str
        Compact label: '1' through '5', and optionally 'TS' or 'EX'
    """
    if not category_str:
        return ''

    category_lower = category_str.lower().strip()

    if category_lower.startswith('cat'):
        try:
            cat_num = int(category_lower[3:])
            if 1 <= cat_num <= 5:
                return str(cat_num)
        except ValueError:
            pass

    match = re.search(r'cat(\d)', category_lower)
    if match:
        cat_num = int(match.group(1))
        if 1 <= cat_num <= 5:
            return str(cat_num)

    if not include_non_hurricane:
        return ''

    if category_lower in ['ts', 'ex']:
        return category_str.upper()

    if 'ts' in category_lower:
        return 'TS'
    if 'ex' in category_lower:
        return 'EX'

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


def plot_max_ele_2d(
        fig, ax,
        maxele_file,
        variable='zeta_max',
        vmin=None, vmax=None, cmap='viridis',
        drawmesh=False, levels=20,
        title=None,
        xmin=None, xmax=None, ymin=None, ymax=None,
        cbar_label=None, cbar_increment=None,
        track_file=None, track_color='red', track_linewidth=2.0, track_markersize=5.0,
        track_annotate_datetime=False, track_annotate_category=False,
        track_annotate_category_inside_circle=False,
        track_annotate_non_hurricane_inside_circle=False, track_annotate_fontsize=8.0,
        draw_shorelines=False, cbar_ticks_increment=None,
        compute_departure=False, departure_reference_file=None, departure_reference_time=None,
        departure_reference_variable=None):
    """
    Plot 2D maximum water level fields from CG ADCIRC maxele NetCDF files.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to plot on
    ax : matplotlib.axes.Axes
        Axes object to plot on
    maxele_file : str
        Path to the maxele NetCDF file (e.g., maxele.63.nc)
    variable : str, optional
        Variable to plot (default: 'zeta_max')
    vmin : float, optional
        Minimum value for colorbar (auto if not specified)
    vmax : float, optional
        Maximum value for colorbar (auto if not specified)
    cmap : str, optional
        Colormap to use (default: 'viridis')
    drawmesh : bool, optional
        Draw mesh lines (default: False)
    levels : int, optional
        Number of contour levels (default: 20)
    title : str, optional
        Custom title for the plot
    xmin : float, optional
        Minimum x value (longitude) for plot range
    xmax : float, optional
        Maximum x value (longitude) for plot range
    ymin : float, optional
        Minimum y value (latitude) for plot range
    ymax : float, optional
        Maximum y value (latitude) for plot range
    cbar_label : str, optional
        Custom label for colorbar (default: auto-generated based on variable)
    cbar_increment : float, optional
        Increment for colorbar levels (overrides levels parameter if specified)
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
    track_annotate_category_inside_circle : bool, optional
        Annotate track points with hurricane category numbers centered inside markers (default: False)
    track_annotate_non_hurricane_inside_circle : bool, optional
        Annotate track points with non-hurricane categories (TS, EX) inside markers (default: False)
    track_annotate_fontsize : float, optional
        Font size for track annotation text (default: 8.0)
    draw_shorelines : bool, optional
        If True, draw 0 m depth contour lines in gray (default: False)
    cbar_ticks_increment : float, optional
        Increment for colorbar tick labels (e.g., 0.5 creates ticks at vmin, vmin+0.5, ..., vmax)
    compute_departure : bool, optional
        If True, compute departure field by subtracting reference water level (default: False)
    departure_reference_file : str, optional
        Path to reference solution file (e.g., fort.63.nc) for departure computation
    departure_reference_time : str, optional
        Reference time in format "YYYY-MM-DD HH:mm:ss" (required if compute_departure=True)
    departure_reference_variable : str, optional
        Variable name to use from the reference file (default: 'zeta')
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    
    # Check if files exist
    if not os.path.exists(maxele_file):
        print(f"Error: File {maxele_file} not found.")
        return False
    
    # Read the mesh data
    print(f"Reading maxele data from {maxele_file}")
    try:
        ds = xr.open_dataset(maxele_file)
    except Exception as e:
        print(f"Error reading {maxele_file}: {e}")
        return False
    
    # Extract mesh information
    adc_x = ds.x.values
    adc_y = ds.y.values
    adc_e = ds.element.values - 1  # Convert to 0-based indexing
    
    # Check if the specified variable exists
    if variable not in ds.variables:
        print(f"Error: Variable '{variable}' not found in {maxele_file}")
        print(f"Available variables: {list(ds.variables.keys())}")
        ds.close()
        return False
    
    # Extract the variable data (maxele data has no time dimension)
    var_shape = ds[variable].shape
    if len(var_shape) == 1:
        # Variables with shape (node,) - defined on nodes
        var_data = ds[variable].values
    elif len(var_shape) == 2 and var_shape[1] == 1:
        # Variables with shape (element, 1) - defined on elements
        var_data = ds[variable].values[:, 0]
    else:
        # For other cases, try to get the data directly
        var_data = ds[variable].values
    
    # Read depth data if needed for shorelines or departure
    depth = None
    if draw_shorelines or compute_departure:
        # Check if depth variable exists
        if 'depth' not in ds.variables:
            if compute_departure:
                print(f"Error: 'depth' variable not found in {maxele_file} (required for departure computation)")
                ds.close()
                return False
            elif draw_shorelines:
                print(f"Error: 'depth' variable not found in {maxele_file} (required for shoreline drawing)")
                ds.close()
                return False
        
        # Extract depth data (should be on nodes, no time dimension)
        depth = ds['depth'].values
    
    # Compute departure if requested
    if compute_departure:
        if departure_reference_file is None or departure_reference_time is None:
            print(f"Error: departure_reference_file and departure_reference_time must be specified when compute_departure=True")
            ds.close()
            return False
        
        if not os.path.exists(departure_reference_file):
            print(f"Error: Reference file {departure_reference_file} not found.")
            ds.close()
            return False
        
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
    
    # Close the dataset
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
        # Mask triangles where any vertex has non-finite value
        mask = np.any(~np.isfinite(var_data[adc_e]), axis=1)
    else:  # Elemental data
        # For elemental data, mask elements with non-finite values
        mask = ~np.isfinite(var_data)
    
    triang.set_mask(mask)
    
    # Create contour plot
    if cbar_increment is not None:
        # Use specified increment for levels
        levels_plot = np.arange(vmin_plot, vmax_plot + cbar_increment, cbar_increment)
        if len(levels_plot) > 0 and abs(levels_plot[-1] - vmax_plot) > cbar_increment * 0.01:
            levels_plot = np.append(levels_plot, vmax_plot)
        levels_plot = levels_plot[levels_plot <= vmax_plot + 1e-9]
        if len(levels_plot) == 0 or abs(levels_plot[-1] - vmax_plot) > 1e-9:
            levels_plot = np.append(levels_plot, vmax_plot)
    else:
        # Use specified number of levels
        levels_plot = np.linspace(vmin_plot, vmax_plot, levels)
    
    contour = ax.tricontourf(triang, var_data, levels=levels_plot, cmap=cmap, extend='both')
    
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
            ax.tricontour(depth_triang, depth, levels=[0.0], colors='black', linewidths=0.2, linestyles='-', zorder=5)
    
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
                if (track_annotate_datetime or track_annotate_category or
                        track_annotate_category_inside_circle):
                    inside_circle_fontsize = min(track_annotate_fontsize, track_markersize * 0.85)
                    for lon, lat, dt, cat in zip(track_lons, track_lats, track_dts, track_categories):
                        if track_annotate_category_inside_circle:
                            category_number = get_category_number(
                                cat, include_non_hurricane=track_annotate_non_hurricane_inside_circle)
                            if category_number:
                                ax.text(lon, lat, category_number, ha='center', va='center',
                                        fontsize=inside_circle_fontsize, color='white',
                                        fontweight='bold', zorder=12)

                        annotation_parts = []

                        # Add datetime if requested
                        if track_annotate_datetime and dt:
                            # Format datetime to short format
                            dt_formatted = format_datetime_short(dt)
                            annotation_parts.append(dt_formatted)

                        # Add category as offset label if requested (not inside-circle mode)
                        if track_annotate_category and not track_annotate_category_inside_circle:
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
    elif variable == 'zeta_max':
        cbar.ax.set_ylabel('Maximum Water Level (m)', rotation=270, labelpad=15)
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
    else:
        ticks = list(cbar.get_ticks())

    # Matplotlib's default tick locator omits vmax; always include colorbar endpoints
    for endpoint in [vmin_plot, vmax_plot]:
        if not any(np.isclose(ticks, endpoint, rtol=0, atol=max(abs(endpoint) * 1e-6, 1e-9))):
            ticks = np.append(ticks, endpoint)
    cbar.set_ticks(np.sort(ticks))
    
    # Set labels and title
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title(f'{variable} from {os.path.basename(maxele_file)}')
    
    # Set plot limits if specified
    if xmin is not None or xmax is not None:
        xlim = [xmin if xmin is not None else ax.get_xlim()[0],
                xmax if xmax is not None else ax.get_xlim()[1]]
        ax.set_xlim(xlim)
    
    if ymin is not None or ymax is not None:
        ylim = [ymin if ymin is not None else ax.get_ylim()[0],
                ymax if ymax is not None else ax.get_ylim()[1]]
        ax.set_ylim(ylim)
    
    # Set aspect ratio
    ax.set_aspect('equal')
    
    return True


def get_parser(add_help=True):
    """Get command line argument parser."""
    parser = argparse.ArgumentParser(description='Plot CG ADCIRC maximum water level fields from maxele NetCDF files.', add_help=add_help)
    parser.add_argument('--maxele-file', type=str, required=True, help='Path to maxele NetCDF file (e.g., maxele.63.nc)')
    parser.add_argument('--variable', type=str, default='zeta_max', help='Variable to plot (default: zeta_max)')
    parser.add_argument('--output', type=str, required=True, help='Output PNG file name')
    parser.add_argument('--figsizex', type=float, default=12.0, help='Figure size in x direction')
    parser.add_argument('--figsizey', type=float, default=10.0, help='Figure size in y direction')
    parser.add_argument('--vmin', type=float, help='Minimum value of the colorbar (auto if not specified)')
    parser.add_argument('--vmax', type=float, help='Maximum value of the colorbar (auto if not specified)')
    parser.add_argument('--cmap', type=str, default='viridis', help='Colormap to use (default: viridis)')
    parser.add_argument('--drawmesh', action='store_true', help='Draw mesh lines')
    parser.add_argument('--levels', type=int, default=20, help='Number of contour levels (default: 20)')
    parser.add_argument('--title', type=str, help='Custom title for the plot')
    parser.add_argument('--xmin', type=float, help='Minimum x value (longitude) for plot range')
    parser.add_argument('--xmax', type=float, help='Maximum x value (longitude) for plot range')
    parser.add_argument('--ymin', type=float, help='Minimum y value (latitude) for plot range')
    parser.add_argument('--ymax', type=float, help='Maximum y value (latitude) for plot range')
    parser.add_argument('--cbar-label', type=str, help='Custom label for colorbar (default: auto-generated based on variable)')
    parser.add_argument('--cbar-increment', type=float, help='Increment for colorbar levels (overrides --levels if specified)')
    parser.add_argument('--cbar-ticks-increment', type=float, help='Increment for colorbar tick labels (e.g., 0.5 creates ticks at vmin, vmin+0.5, ..., vmax)')
    parser.add_argument('--track-file', type=str, help='Path to KMZ file containing hurricane best track to overlay on plot')
    parser.add_argument('--track-color', type=str, default='red', help='Color for the track line (default: red)')
    parser.add_argument('--track-linewidth', type=float, default=2.0, help='Line width for the track (default: 2.0)')
    parser.add_argument('--track-markersize', type=float, default=5.0, help='Marker size for track points (default: 5.0)')
    parser.add_argument('--track-annotate-datetime', action='store_true', help='Annotate track points with datetime labels')
    parser.add_argument('--track-annotate-category', action='store_true', help='Annotate track points with category labels')
    parser.add_argument('--track-annotate-category-inside-circle', action='store_true', help='Annotate track points with hurricane category numbers centered inside markers')
    parser.add_argument('--track-annotate-non-hurricane-inside-circle', action='store_true', help='Annotate track points with non-hurricane categories (TS, EX) inside markers (requires --track-annotate-category-inside-circle)')
    parser.add_argument('--track-annotate-fontsize', type=float, default=8.0, help='Font size for track annotation text (default: 8.0)')
    parser.add_argument('--draw-shorelines', action='store_true', help='Draw 0 m depth contour lines in gray')
    parser.add_argument('--departure', action='store_true', help='Compute and plot departure field by subtracting reference water level')
    parser.add_argument('--departure-reference-file', type=str, help='Path to reference solution file (e.g., fort.63.nc) for departure computation')
    parser.add_argument('--departure-reference-time', type=str, help='Reference time in format "YYYY-MM-DD HH:mm:ss" (required if --departure is used)')
    parser.add_argument('--departure-reference-variable', type=str, default='zeta', help='Variable name to use from the reference file (default: zeta)')
    
    return parser


def main(args=None):
    """Main function for command line usage."""
    if args is None:
        args = get_parser().parse_args()
    
    # Validate departure options
    if args.departure:
        if args.departure_reference_file is None or args.departure_reference_time is None:
            print("Error: --departure-reference-file and --departure-reference-time must be specified when using --departure")
            sys.exit(1)
    
    # Create the plot
    print(f"Creating maxele plot")
    fig, ax = plt.subplots(figsize=(args.figsizex, args.figsizey))
    
    success = plot_max_ele_2d(
        fig, ax,
        maxele_file=args.maxele_file,
        variable=args.variable,
        vmin=args.vmin, vmax=args.vmax, cmap=args.cmap,
        drawmesh=args.drawmesh, levels=args.levels,
        title=args.title,
        xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax,
        cbar_label=args.cbar_label, cbar_increment=args.cbar_increment,
        track_file=args.track_file, track_color=args.track_color,
        track_linewidth=args.track_linewidth, track_markersize=args.track_markersize,
        track_annotate_datetime=args.track_annotate_datetime,
        track_annotate_category=args.track_annotate_category,
        track_annotate_category_inside_circle=args.track_annotate_category_inside_circle,
        track_annotate_non_hurricane_inside_circle=args.track_annotate_non_hurricane_inside_circle,
        track_annotate_fontsize=args.track_annotate_fontsize,
        draw_shorelines=args.draw_shorelines, cbar_ticks_increment=args.cbar_ticks_increment,
        compute_departure=args.departure, departure_reference_file=args.departure_reference_file,
        departure_reference_time=args.departure_reference_time,
        departure_reference_variable=args.departure_reference_variable
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
