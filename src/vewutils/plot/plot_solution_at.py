def plot_solution_at(
        fig, ax,
        lon, lat,
        date_start, date_end, 
        solution_files, solution_starts, solution_labels, 
        var_names, ylabel, title,
        plot_movingaverage=False, solution_colors=None, plot_in_foot=False, 
        movingaverage_window=0, ymin=None, ymax=None, mag=False, options=None):
    """
    Plot solution variables at a specific point from NetCDF solution files.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to plot on
    ax : matplotlib.axes.Axes
        Axes object to plot on
    lon : float
        Longitude of the point
    lat : float
        Latitude of the point
    date_start : datetime
        Start date for plotting
    date_end : datetime
        End date for plotting
    solution_files : list
        List of solution file paths or list of lists for multiple files per solution
    solution_starts : list
        List of start times for each solution file
    solution_labels : list
        List of labels for each solution
    var_names : str or list of str
        Variable name(s) to plot (e.g., 'zeta', ['u-vel', 'v-vel'], etc.)
    ylabel : str
        Label for y-axis
    title : str
        Title for the plot
    plot_movingaverage : bool, optional
        Whether to plot moving average, by default False
    solution_colors : str or list, optional
        Color(s) for solutions. If single color string, used for all solutions.
        If list, colors for each solution. By default None (auto-colors)
    plot_in_foot : bool, optional
        Whether to convert units to feet, by default False
    movingaverage_window : int, optional
        Window size for moving average, by default 0
    ymin : float, optional
        Minimum y-axis value, by default None (auto)
    ymax : float, optional
        Maximum y-axis value, by default None (auto)
    mag : bool, optional
        Compute and plot magnitude from two variables, by default False
    options : dict, optional
        Additional options, by default None
    """
    import os
    import sys
    import requests
    from datetime import datetime, timedelta
    import matplotlib.pyplot as plt
    import pytz
    from matplotlib.dates import DateFormatter
    import pandas as pd
    import numpy as np
    from vewutils.plot.get_solution_at import get_solution_at
    
    # Handle single variable name as string
    if isinstance(var_names, str):
        var_names = [var_names]
    
    # Set unit conversion factor
    if plot_in_foot:
        m2ft = 1.0/0.3048
    else:
        m2ft = 1.0
        
    # Get the solution data
    solution_times = []
    solution_data = {var: [] for var in var_names}
    
    for i in range(len(solution_files)):
        if isinstance(solution_files[i], list):
            solution_file = solution_files[i]
            solution_start = solution_starts[i]
        else:
            solution_file = [solution_files[i]]
            solution_start = [solution_starts[i]]
        
        solution_time = []
        solution_vars = {var: [] for var in var_names}
        
        print('Processing solution files:')
        for j in range(len(solution_file)):
            solution_filej = solution_file[j]
            solution_startj = solution_start[j]
            print(f'  Processing: {solution_filej}')
            
            # Get solution data using get_solution_at
            result = get_solution_at(solution_filej, lon, lat, var_names)
            solution_timej = result['time']
            
            if solution_startj:
                tdj = solution_startj - solution_timej[0]
                solution_timej = [tdj + t for t in solution_timej]
            
            if solution_timej is not None:
                solution_time.extend(solution_timej.tolist())
                for var in var_names:
                    if var in result:
                        solution_vars[var].extend(result[var].tolist())
                    else:
                        solution_vars[var].extend([np.nan] * len(solution_timej))
        
        solution_times.append(solution_time)
        for var in var_names:
            # Filter out extreme values
            var_data = [val if abs(val) <= 100 else np.nan for val in solution_vars[var]]
            solution_data[var].append(var_data)
        print('  Done.')
    
    print(f'Processed {len(solution_files)} solution file(s) successfully.')
    
    # Compute magnitude if requested
    if mag:
        if len(var_names) != 2:
            raise ValueError("Magnitude computation requires exactly 2 variables")
        
        print('Computing magnitude from two variables...')
        # Create magnitude data for each solution
        for i in range(len(solution_files)):
            var1_data = solution_data[var_names[0]][i]
            var2_data = solution_data[var_names[1]][i]
            
            # Compute magnitude: sqrt(var1^2 + var2^2)
            magnitude_data = []
            for j in range(len(var1_data)):
                if (not np.isnan(var1_data[j]) and not np.isnan(var2_data[j]) and 
                    var1_data[j] is not None and var2_data[j] is not None):
                    mag_val = np.sqrt(var1_data[j]**2 + var2_data[j]**2)
                    magnitude_data.append(mag_val)
                else:
                    magnitude_data.append(np.nan)
            
            # Add magnitude to solution data
            solution_data['magnitude'] = solution_data.get('magnitude', [])
            solution_data['magnitude'].append(magnitude_data)
        
        # Update var_names to include magnitude
        var_names = ['magnitude']
        print('  Done.')
    
    # Calculate moving average over a 2-day window
    if plot_movingaverage:
        window_size = timedelta(days=2)
        
        solution_times_ma = []
        solution_data_ma = {var: [] for var in var_names}
        
        for i in range(len(solution_files)):
            times = solution_times[i]
            times_ma = []
            vars_ma = {var: [] for var in var_names}
            
            for j in range(len(times)):
                start_time = times[j] - window_size / 2
                end_time = times[j] + window_size / 2
                
                start_i = 0
                end_i = 0
                for start_i in reversed(range(j)):
                    if times[start_i] < start_time:
                        start_i = min(j, start_i + 1)
                        break
                for end_i in range(j, len(times)):
                    if times[end_i] > end_time:
                        end_i = max(j, end_i - 1)
                        break
                
                for var in var_names:
                    var_window = [solution_data[var][i][k] for k in range(start_i, end_i+1)]
                    var_window_valid = [v for v in var_window if v is not None and not (isinstance(v, float) and np.isnan(v))]
                    if var_window_valid:
                        var_avg = sum(var_window_valid) / len(var_window_valid)
                    else:
                        var_avg = np.nan
                    vars_ma[var].append(var_avg)
                
                times_ma.append(times[j])
            
            solution_times_ma.append(times_ma)
            for var in var_names:
                solution_data_ma[var].append(vars_ma[var])
    
    # Plot the solution data for each variable
    for var_idx, var_name in enumerate(var_names):
        for i in range(len(solution_files)):
            if solution_colors is not None:
                if isinstance(solution_colors, str):
                    # Single color for all solutions
                    solution_color = solution_colors
                else:
                    # List of colors for each solution
                    solution_color = solution_colors[i % len(solution_colors)]
            else:
                # Cycle through colors
                colors = ['b', 'r', 'y', 'g', 'c', 'm', 'k']
                solution_color = colors[i % len(colors)]
            
            # Apply moving average if requested
            var_data = solution_data[var_name][i]
            if movingaverage_window > 0:
                var_data_ma = []
                for j in range(len(var_data)):
                    start = max(0, j-1)
                    end = min(len(var_data)-1, j+1)
                    window = var_data[start:end+1]
                    window_valid = [v for v in window if v is not None and not (isinstance(v, float) and np.isnan(v))]
                    if window_valid:
                        ma_val = np.nanmean(window_valid)
                    else:
                        ma_val = np.nan
                    var_data_ma.append(ma_val)
            else:
                var_data_ma = var_data
            
            # Create label with variable name if multiple variables
            if len(var_names) > 1:
                label = f'{solution_labels[i]} ({var_name})'
            else:
                label = solution_labels[i]
            
            ax.plot(solution_times[i], np.array(var_data_ma) * m2ft, '-', 
                   color=solution_color, label=label, alpha=0.8)
    
    # Plot moving averages if requested
    if plot_movingaverage:
        fmt = '--'
        for var_idx, var_name in enumerate(var_names):
            for i in range(len(solution_files)):
                if solution_colors is not None:
                    if isinstance(solution_colors, str):
                        # Single color for all solutions
                        solution_color = solution_colors
                    else:
                        # List of colors for each solution
                        solution_color = solution_colors[i % len(solution_colors)]
                else:
                    colors = ['b', 'r', 'y', 'g', 'c', 'm', 'k']
                    solution_color = colors[i % len(colors)]
                
                if len(var_names) > 1:
                    label = f'{solution_labels[i]} ({var_name}) (2d MA)'
                else:
                    label = f'{solution_labels[i]} (2d MA)'
                
                ax.plot(solution_times_ma[i], np.array(solution_data_ma[var_name][i]) * m2ft, 
                       fmt, color=solution_color, label=label, alpha=0.6)
    
    # Format the plot
    ax.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M'))
    fig.autofmt_xdate()
    ax.grid(True, alpha=0.3)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlim([date_start, date_end])
    
    # Set y-axis limits if specified
    if ymin is not None or ymax is not None:
        ax.set_ylim(bottom=ymin, top=ymax)
    
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    return True

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(add_help=False, description='Plot solution variables at a point.')
    parser.add_argument('--lon', type=float, required=True, help='Longitude')
    parser.add_argument('--lat', type=float, required=True, help='Latitude')
    parser.add_argument('--date-start', type=str, required=True, help='Start date (YYYY-MM-DD)')
    parser.add_argument('--date-end', type=str, required=True, help='End date (YYYY-MM-DD)')
    parser.add_argument('--solution-files', type=str, nargs='+', required=True, help='List of solution files')
    parser.add_argument('--solution-starts', type=str, nargs='*', help='List of solution start times (default: None for all files)')
    parser.add_argument('--solution-labels', type=str, nargs='*', help='List of solution labels (default: auto-generated)')
    parser.add_argument('--var-names', type=str, nargs='+', required=True, help='Variable name(s) to plot')
    parser.add_argument('--ylabel', type=str, help='Y-axis label (default: auto-generated based on variables)')
    parser.add_argument('--title', type=str, help='Plot title (default: auto-generated)')
    parser.add_argument('--plot-movingaverage', action='store_true', help='Plot moving average')
    parser.add_argument('--plot-in-foot', action='store_true', help='Convert units to feet')
    parser.add_argument('--movingaverage-window', type=int, default=0, help='Moving average window size')
    parser.add_argument('--solution-colors', type=str, nargs='*', help='Color(s) for solutions. Single color for all, or list of colors for each solution')
    parser.add_argument('--fig-width', type=float, default=10.0, help='Figure width in inches (default: 10.0)')
    parser.add_argument('--fig-height', type=float, default=6.0, help='Figure height in inches (default: 6.0)')
    parser.add_argument('--ymin', type=float, help='Minimum y-axis value (default: auto)')
    parser.add_argument('--ymax', type=float, help='Maximum y-axis value (default: auto)')
    parser.add_argument('--mag', action='store_true', help='Compute and plot magnitude from two variables (requires exactly 2 var-names)')
    parser.add_argument('--outputfile', type=str, required=True, help='Output figure file name')
    return parser

def main(args=None):
    from datetime import datetime
    import matplotlib.pyplot as plt
    if args is None:
        args = get_parser().parse_args()
    
    # Validate --mag option requires exactly 2 variables
    if args.mag and len(args.var_names) != 2:
        raise ValueError("--mag option requires exactly 2 variable names")
    
    date_start = datetime.strptime(args.date_start, '%Y-%m-%d')
    date_end = datetime.strptime(args.date_end, '%Y-%m-%d')
    
    # Handle optional solution_starts
    if args.solution_starts:
        solution_starts = [datetime.strptime(sol_start, '%Y-%m-%d') if sol_start else None 
                          for sol_start in args.solution_starts]
    else:
        solution_starts = [None] * len(args.solution_files)
    
    # Handle optional solution_labels
    if args.solution_labels:
        solution_labels = args.solution_labels
    else:
        solution_labels = [f"Solution {i+1}" for i in range(len(args.solution_files))]
    
    # Handle optional ylabel
    if args.ylabel:
        ylabel = args.ylabel
    else:
        # Auto-generate ylabel based on variables
        if args.mag:
            ylabel = 'Magnitude (units)'
        elif len(args.var_names) == 1:
            var_name = args.var_names[0]
            if var_name == 'zeta':
                ylabel = 'Water Level (m)'
            elif var_name in ['u-vel', 'v-vel']:
                ylabel = 'Velocity (m/s)'
            else:
                ylabel = f'{var_name} (units)'
        else:
            ylabel = 'Variable Values (units)'
    
    # Handle optional title
    if args.title:
        title = args.title
    else:
        # Auto-generate title
        if args.mag:
            var_str = f'Magnitude of {args.var_names[0]} and {args.var_names[1]}'
        else:
            var_str = ', '.join(args.var_names)
        title = f'{var_str} at Point ({args.lon}, {args.lat})'
    
    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    plot_solution_at(
        fig, ax,
        args.lon, args.lat,
        date_start, date_end, 
        args.solution_files, solution_starts, solution_labels,
        args.var_names, ylabel, title,
        args.plot_movingaverage, solution_colors=args.solution_colors, 
        plot_in_foot=args.plot_in_foot, movingaverage_window=args.movingaverage_window,
        ymin=args.ymin, ymax=args.ymax, mag=args.mag)
    
    plt.tight_layout()
    plt.savefig(args.outputfile, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {args.outputfile}")

if __name__ == '__main__':
    main()
