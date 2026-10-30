_LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]


def _get_solution_color(index, solution_colors=None):
    if solution_colors is not None:
        if isinstance(solution_colors, str):
            return solution_colors
        return solution_colors[index % len(solution_colors)]

    colors = ["b", "r", "y", "g", "c", "m", "k"]
    return colors[index % len(colors)]


def _get_solution_linestyle(index, distinct_linestyles=True, n_series=1):
    if not distinct_linestyles or n_series <= 1:
        return "-"
    return _LINESTYLES[index % len(_LINESTYLES)]


def _get_unit_scale(plot_in_foot=False, unit_scale=None):
    if unit_scale is not None:
        return unit_scale
    if plot_in_foot:
        return 1.0 / 0.3048
    return 1.0


def plot_solution_along_transect(
    fig,
    ax,
    distance,
    series,
    ylabel,
    title,
    plot_in_foot=False,
    unit_scale=None,
    ymin=None,
    ymax=None,
    solution_colors=None,
    distinct_linestyles=True,
    secondary_values=None,
    secondary_label=None,
    secondary_color=None,
):
    """
    Plot solution variable profile(s) along a transect.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to plot on
    ax : matplotlib.axes.Axes
        Axes object to plot on
    distance : array-like
        Distance along transect in meters
    series : list of dict
        Each item contains ``values`` and ``label``, and may contain ``color``
        and ``linestyle``
    ylabel : str
        Label for y-axis
    title : str
        Title for the plot
    plot_in_foot : bool, optional
        Whether to convert units to feet, by default False
    unit_scale : float, optional
        Explicit scale factor applied to plotted values. If None, uses meters-to-feet
        conversion when plot_in_foot is True.
    ymin, ymax : float, optional
        Y-axis limits
    solution_colors : str or list of str, optional
        Color(s) for solution series
    distinct_linestyles : bool, optional
        When True and more than one solution series is plotted, cycle line styles
        so each solution is distinguishable. By default True.
    secondary_values : array-like, optional
        Secondary variable values to plot on the same axes
    secondary_label : str, optional
        Label for the secondary line
    secondary_color : str, optional
        Color for the secondary line
    """
    import numpy as np

    scale = _get_unit_scale(plot_in_foot=plot_in_foot, unit_scale=unit_scale)

    n_series = len(series)
    for i, item in enumerate(series):
        color = item.get("color")
        if color is None:
            color = _get_solution_color(i, solution_colors=solution_colors)
        linestyle = item.get("linestyle")
        if linestyle is None:
            linestyle = _get_solution_linestyle(
                i, distinct_linestyles=distinct_linestyles, n_series=n_series
            )

        ax.plot(
            distance,
            np.array(item["values"]) * scale,
            color=color,
            linestyle=linestyle,
            label=item["label"],
            alpha=0.8,
        )

    if secondary_values is not None:
        ax.plot(
            distance,
            np.array(secondary_values) * scale,
            "-",
            color=secondary_color or "saddlebrown",
            label=secondary_label,
            alpha=0.8,
        )

    ax.grid(True, alpha=0.3)
    ax.set_title(title)
    ax.set_xlabel("Distance along transect (m)")
    ax.set_ylabel(ylabel)

    if ymin is not None or ymax is not None:
        ax.set_ylim(bottom=ymin, top=ymax)

    if series or secondary_label is not None:
        ax.legend(loc="best")

    return True


def _format_time_label(selected_time):
    return selected_time.strftime("%Y-%m-%d %H:%M:%S UTC")


def _transect_title_suffix(start_x, start_y, end_x, end_y):
    return f"Transect ({start_x:g}, {start_y:g}) to ({end_x:g}, {end_y:g})"


PLOT_VARIABLES = ["zeta", "u-vel", "v-vel", "v-tangential", "flux"]


class _FramePlot:
    """
    One figure reused for every frame of a plot: the first frame is drawn with
    plot_solution_along_transect, later frames only update line data and title.
    """

    def __init__(self, fig_width, fig_height, plot_kwargs):
        self.figsize = (fig_width, fig_height)
        self.plot_kwargs = plot_kwargs
        self.fixed_ylim = (
            plot_kwargs.get("ymin") is not None and plot_kwargs.get("ymax") is not None
        )
        self.fig = None
        self.bbox = None

    def save(self, outputfile, distance, series, title, secondary_values=None):
        import numpy as np

        if self.fig is None:
            import matplotlib.pyplot as plt

            self.fig, self.ax = plt.subplots(figsize=self.figsize)
            plot_solution_along_transect(
                self.fig,
                self.ax,
                distance,
                series,
                title=title,
                secondary_values=secondary_values,
                **self.plot_kwargs,
            )
            self.lines = list(self.ax.get_lines())
            self.scale = _get_unit_scale(
                plot_in_foot=self.plot_kwargs.get("plot_in_foot", False),
                unit_scale=self.plot_kwargs.get("unit_scale"),
            )
        else:
            values = [item["values"] for item in series]
            if secondary_values is not None:
                values.append(secondary_values)
            for line, v in zip(self.lines, values):
                line.set_ydata(np.array(v) * self.scale)
            self.ax.set_title(title)
            if self.bbox is not None:
                # fixed axis limits: layout and crop box are those of the first frame
                self.fig.savefig(outputfile, dpi=150, bbox_inches=self.bbox)
                print(f"Plot saved to {outputfile}")
                return
            kw = self.plot_kwargs
            self.ax.relim()
            self.ax.autoscale_view(scalex=False)
            if kw.get("ymin") is not None or kw.get("ymax") is not None:
                self.ax.set_ylim(bottom=kw.get("ymin"), top=kw.get("ymax"))
        self.fig.tight_layout()
        self.fig.savefig(outputfile, dpi=150, bbox_inches="tight")
        if self.fixed_ylim and self.bbox is None:
            import matplotlib as mpl

            # measured at the saving dpi, as savefig does
            dpi = self.fig.dpi
            self.fig.dpi = 150
            renderer = self.fig.canvas.get_renderer()
            self.bbox = self.fig.get_tightbbox(renderer).padded(
                mpl.rcParams["savefig.pad_inches"]
            )
            self.fig.dpi = dpi
        print(f"Plot saved to {outputfile}")

    def close(self):
        if self.fig is not None:
            import matplotlib.pyplot as plt

            plt.close(self.fig)


def _output_prefix(template, timestep, frame, multi):
    if not multi:
        return template
    if "{" in template:
        return template.format(timestep=timestep, frame=frame)
    return f"{template}_{timestep:05d}"


def _resolve_solution_labels(solution_labels, n_solutions):
    if solution_labels:
        if len(solution_labels) != n_solutions:
            raise ValueError(
                f"Number of --solution-labels ({len(solution_labels)}) must match "
                f"number of solutions ({n_solutions})"
            )
        return solution_labels
    return [f"Solution {i + 1}" for i in range(n_solutions)]


def _extract_transect_solutions(
    solution_files,
    start_x,
    start_y,
    end_x,
    end_y,
    spacing,
    var_names,
    timestep,
    datetime_str,
    use_cache,
    tvv_fort142=None,
    tvv_search_radius=None,
):
    from vewutils.plot.get_solution_along_transect import get_solution_along_transect

    solutions = []
    for i, solution_file in enumerate(solution_files):
        print(f"  Processing: {solution_file}")
        solutions.append(
            get_solution_along_transect(
                solution_file,
                start_x,
                start_y,
                end_x,
                end_y,
                spacing,
                var_names,
                timestep=timestep,
                datetime_str=datetime_str,
                use_cache=use_cache,
                # only the first solution's depth is plotted, so only it is moved
                tvv_fort142=tvv_fort142 if i == 0 else None,
                tvv_search_radius=tvv_search_radius,
            )
        )
    return solutions


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        add_help=False,
        description="Plot solution profiles along a transect at a single time step.",
    )
    parser.add_argument("--start-x", type=float, required=True, help="Start x-coordinate")
    parser.add_argument("--start-y", type=float, required=True, help="Start y-coordinate")
    parser.add_argument("--end-x", type=float, required=True, help="End x-coordinate")
    parser.add_argument("--end-y", type=float, required=True, help="End y-coordinate")
    parser.add_argument(
        "--spacing",
        type=float,
        required=True,
        help="Spacing between extraction points in meters",
    )
    time_group = parser.add_mutually_exclusive_group(required=True)
    time_group.add_argument(
        "--timestep",
        type=int,
        help="Time step to plot (0-based index, use -1 for last time step)",
    )
    time_group.add_argument(
        "--datetime",
        type=str,
        help="Target time in YYYY-MM-DDThh:mm:ss format",
    )
    time_group.add_argument(
        "--timesteps",
        type=str,
        nargs="+",
        help=(
            "Several time steps to plot, one set of figures each: indices, inclusive "
            "ranges first:last[:step], or comma-separated lists of these "
            "(e.g. 159:1359). Each file is opened once and the transect "
            "interpolation computed once. --output-prefix may contain {timestep} "
            "and {frame} (position in the list plus --frame-offset), e.g. "
            "frames/zeta_{frame:04d}; otherwise _{timestep:05d} is appended."
        ),
    )
    parser.add_argument(
        "--fort63",
        type=str,
        nargs="+",
        help="Path(s) to fort.63.nc file(s) for water level extraction",
    )
    parser.add_argument(
        "--fort64",
        type=str,
        nargs="+",
        help="Path(s) to fort.64.nc file(s) for velocity extraction",
    )
    parser.add_argument(
        "--solution-labels",
        type=str,
        nargs="*",
        help="List of solution labels (default: auto-generated)",
    )
    parser.add_argument(
        "--solution-colors",
        type=str,
        nargs="*",
        help="Color(s) for solutions. Single color for all, or list of colors for each solution",
    )
    parser.add_argument(
        "--distinct-linestyles",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Use a different line style for each solution when more than one "
            "solution file is plotted (default: enabled). "
            "Pass --no-distinct-linestyles to draw every solution as a solid line."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        required=True,
        help="Prefix for output figure file names",
    )
    parser.add_argument(
        "--plot-in-foot",
        action="store_true",
        help="Convert elevation and velocity units to feet",
    )
    parser.add_argument(
        "--fig-width",
        type=float,
        default=10.0,
        help="Figure width in inches (default: 10.0)",
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=4.0,
        help="Figure height in inches (default: 6.0)",
    )
    parser.add_argument("--ymin", type=float, help="Minimum y-axis value (default: auto)")
    parser.add_argument("--ymax", type=float, help="Maximum y-axis value (default: auto)")
    parser.add_argument(
        "--fort142",
        type=str,
        help=(
            "Time-varying weir file (fort.142). Its VaryType=4 lines (time-varying "
            "crest VEWs) move the wall-top nodes in the bathymetry line to their "
            "crest at the plotted time. Crest tables are read relative to its directory. "
            "Applies to the first --fort63 file, whose depth is the one plotted."
        ),
    )
    parser.add_argument(
        "--tvv-search-radius",
        type=float,
        help=(
            "Search radius (m) for VaryType=4 lines without SearchRadius, as "
            "TVV_SEARCH_RADIUS in fort.15 (default: 1e-6)"
        ),
    )
    parser.add_argument(
        "--use-cache",
        action="store_true",
        help=(
            "Save extracted transect data to on-disk cache files (default: do not "
            "cache). Not used with --timesteps."
        ),
    )
    parser.add_argument(
        "--variables",
        type=str,
        nargs="+",
        choices=PLOT_VARIABLES,
        help=(
            "Plots to make (default: every plot the given files allow). "
            "zeta needs --fort63, the velocities --fort64, flux both."
        ),
    )
    parser.add_argument(
        "--frame-offset",
        type=int,
        default=0,
        help="Added to the list position for {frame} in --output-prefix (default: 0)",
    )
    return parser


def main(args=None):
    import numpy as np

    if args is None:
        args = get_parser().parse_args()

    fort63_files = args.fort63 or []
    fort64_files = args.fort64 or []

    if not fort63_files and not fort64_files:
        raise ValueError("At least one of --fort63 or --fort64 must be specified")

    if args.spacing <= 0:
        raise ValueError("--spacing must be positive")

    if fort63_files and fort64_files and len(fort63_files) != len(fort64_files):
        raise ValueError(
            "When both --fort63 and --fort64 are specified, the number of files must match"
        )

    variables = getattr(args, "variables", None)
    if variables is None:
        variables = []
        if fort63_files:
            variables.append("zeta")
        if fort64_files:
            variables += ["u-vel", "v-vel", "v-tangential"]
        if fort63_files and fort64_files:
            variables.append("flux")
    for var in variables:
        if var in ("zeta", "flux") and not fort63_files:
            raise ValueError(f"--variables {var} needs --fort63")
        if var != "zeta" and not fort64_files:
            raise ValueError(f"--variables {var} needs --fort64")
    need63 = "zeta" in variables or "flux" in variables
    need64 = any(var != "zeta" for var in variables)

    n_solutions = max(len(fort63_files), len(fort64_files))
    solution_labels = _resolve_solution_labels(args.solution_labels, n_solutions)
    solution_colors = args.solution_colors
    if solution_colors is not None and len(solution_colors) not in (1, n_solutions):
        raise ValueError(
            "Provide either one --solution-colors value for all solutions or one per solution"
        )

    transect_suffix = _transect_title_suffix(
        args.start_x, args.start_y, args.end_x, args.end_y
    )
    timesteps_arg = getattr(args, "timesteps", None)
    multi = timesteps_arg is not None
    transect = (args.start_x, args.start_y, args.end_x, args.end_y, args.spacing)

    # zeta_frames[k][i]: solution i at frame k (same for velocity_frames)
    zeta_frames = []
    velocity_frames = []

    if multi:
        from vewutils.plot.get_solution_along_transect import (
            get_solution_along_transect_timesteps,
            parse_timesteps,
        )

        timesteps = parse_timesteps(timesteps_arg)

        def extract(files, var_names, tvv_fort142=None):
            per_file = []
            for i, solution_file in enumerate(files):
                print(f"  Processing: {solution_file}")
                per_file.append(
                    get_solution_along_transect_timesteps(
                        solution_file,
                        *transect,
                        var_names,
                        timesteps,
                        # only the first solution's depth is plotted, so only it is moved
                        tvv_fort142=tvv_fort142 if i == 0 else None,
                        tvv_search_radius=args.tvv_search_radius,
                    )
                )
            return [list(frame) for frame in zip(*per_file)]
    else:
        timesteps = [args.timestep]

        def extract(files, var_names, tvv_fort142=None):
            return [
                _extract_transect_solutions(
                    files,
                    *transect,
                    var_names,
                    args.timestep,
                    args.datetime,
                    args.use_cache,
                    tvv_fort142=tvv_fort142,
                    tvv_search_radius=args.tvv_search_radius,
                )
            ]

    if need63:
        print("Extracting fort.63 data:")
        zeta_frames = extract(fort63_files, ["zeta", "depth"], tvv_fort142=args.fort142)

    if need64:
        print("Extracting fort.64 data:")
        velocity_frames = extract(fort64_files, ["u-vel", "v-vel"])

    dx = args.end_x - args.start_x
    dy = args.end_y - args.start_y
    transect_length = np.hypot(dx, dy)
    if transect_length > 0:
        tx = dx / transect_length
        ty = dy / transect_length
    else:
        tx = 1.0
        ty = 0.0

    elev_unit = "ft" if args.plot_in_foot else "m"
    vel_unit = "ft/s" if args.plot_in_foot else "m/s"
    flux_unit = f"{elev_unit}^2/s"

    common = dict(
        ymin=args.ymin,
        ymax=args.ymax,
        solution_colors=solution_colors,
        distinct_linestyles=args.distinct_linestyles,
    )
    plots = {}
    if "zeta" in variables:
        plots["zeta"] = _FramePlot(
            args.fig_width,
            args.fig_height,
            dict(
                common,
                ylabel=f"Elevation ({elev_unit})",
                plot_in_foot=args.plot_in_foot,
                secondary_label="Bathymetry (-depth)",
                secondary_color="saddlebrown",
            ),
        )
    velocity_plots = [
        ("u-vel", f"U-Velocity ({vel_unit})", "U-Velocity"),
        ("v-vel", f"V-Velocity ({vel_unit})", "V-Velocity"),
        ("v-tangential", f"Tangential Velocity ({vel_unit})", "Tangential Velocity"),
    ]
    for var_key, ylabel, _ in velocity_plots:
        if var_key in variables:
            plots[var_key] = _FramePlot(
                args.fig_width,
                args.fig_height,
                dict(common, ylabel=ylabel, plot_in_foot=args.plot_in_foot),
            )
    if "flux" in variables:
        plots["flux"] = _FramePlot(
            args.fig_width,
            args.fig_height,
            dict(
                common,
                ylabel=f"Flux ({flux_unit})",
                unit_scale=(1.0 / 0.3048) ** 2 if args.plot_in_foot else 1.0,
            ),
        )

    try:
        for k, timestep in enumerate(timesteps):
            zeta_solutions = zeta_frames[k] if zeta_frames else []
            velocity_solutions = velocity_frames[k] if velocity_frames else []
            reference_solution = (
                zeta_solutions[0] if zeta_solutions else velocity_solutions[0]
            )
            time_label = _format_time_label(reference_solution["time"])
            distance = reference_solution["distance"]
            prefix = _output_prefix(
                args.output_prefix,
                reference_solution["timestep"],
                k + args.frame_offset,
                multi,
            )

            if "zeta" in plots:
                zeta_series = [
                    {"values": solution["zeta"], "label": solution_labels[i]}
                    for i, solution in enumerate(zeta_solutions)
                ]
                plots["zeta"].save(
                    f"{prefix}_zeta.png",
                    distance,
                    zeta_series,
                    f"Water Level and Bathymetry at {time_label}\n{transect_suffix}",
                    secondary_values=-np.array(zeta_solutions[0]["depth"]),
                )

            for var_key, _, plot_name in velocity_plots:
                if var_key not in plots:
                    continue
                series = []
                for i, solution in enumerate(velocity_solutions):
                    if var_key == "v-tangential":
                        values = solution["u-vel"] * tx + solution["v-vel"] * ty
                    else:
                        values = solution[var_key]
                    series.append({"values": values, "label": solution_labels[i]})
                plots[var_key].save(
                    f"{prefix}_{var_key}.png",
                    distance,
                    series,
                    f"{plot_name} at {time_label}\n{transect_suffix}",
                )

            if "flux" in plots:
                flux_series = []
                for i, (zeta_solution, velocity_solution) in enumerate(
                    zip(zeta_solutions, velocity_solutions)
                ):
                    v_tangential = (
                        velocity_solution["u-vel"] * tx + velocity_solution["v-vel"] * ty
                    )
                    total_depth = np.array(zeta_solution["zeta"]) + np.array(
                        zeta_solution["depth"]
                    )
                    flux_series.append(
                        {"values": total_depth * v_tangential, "label": solution_labels[i]}
                    )
                plots["flux"].save(
                    f"{prefix}_flux.png",
                    distance,
                    flux_series,
                    f"Flux at {time_label}\n{transect_suffix}",
                )
    finally:
        for plot in plots.values():
            plot.close()


if __name__ == "__main__":
    main()
