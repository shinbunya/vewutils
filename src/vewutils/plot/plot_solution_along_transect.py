def _get_solution_color(index, solution_colors=None):
    if solution_colors is not None:
        if isinstance(solution_colors, str):
            return solution_colors
        return solution_colors[index % len(solution_colors)]

    colors = ["b", "r", "y", "g", "c", "m", "k"]
    return colors[index % len(colors)]


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
    secondary_values : array-like, optional
        Secondary variable values to plot on the same axes
    secondary_label : str, optional
        Label for the secondary line
    secondary_color : str, optional
        Color for the secondary line
    """
    import numpy as np

    scale = _get_unit_scale(plot_in_foot=plot_in_foot, unit_scale=unit_scale)

    for i, item in enumerate(series):
        color = item.get("color")
        if color is None:
            color = _get_solution_color(i, solution_colors=solution_colors)

        ax.plot(
            distance,
            np.array(item["values"]) * scale,
            "-",
            color=color,
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


def _save_figure(fig, outputfile):
    import matplotlib.pyplot as plt

    plt.tight_layout()
    fig.savefig(outputfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {outputfile}")


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
):
    from vewutils.plot.get_solution_along_transect import get_solution_along_transect

    solutions = []
    for solution_file in solution_files:
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
        "--use-cache",
        action="store_true",
        help="Save extracted transect data to on-disk cache files (default: do not cache)",
    )
    return parser


def main(args=None):
    import numpy as np
    import matplotlib.pyplot as plt

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

    n_solutions = max(len(fort63_files), len(fort64_files))
    solution_labels = _resolve_solution_labels(args.solution_labels, n_solutions)
    solution_colors = args.solution_colors
    if solution_colors is not None and len(solution_colors) not in (1, n_solutions):
        raise ValueError(
            "Provide either one --solution-colors value for all solutions or one per solution"
        )

    timestep = args.timestep
    datetime_str = args.datetime
    transect_suffix = _transect_title_suffix(
        args.start_x, args.start_y, args.end_x, args.end_y
    )
    use_cache = args.use_cache

    zeta_solutions = []
    velocity_solutions = []

    if fort63_files:
        print("Extracting fort.63 data:")
        zeta_solutions = _extract_transect_solutions(
            fort63_files,
            args.start_x,
            args.start_y,
            args.end_x,
            args.end_y,
            args.spacing,
            ["zeta", "depth"],
            timestep,
            datetime_str,
            use_cache,
        )

    if fort64_files:
        print("Extracting fort.64 data:")
        velocity_solutions = _extract_transect_solutions(
            fort64_files,
            args.start_x,
            args.start_y,
            args.end_x,
            args.end_y,
            args.spacing,
            ["u-vel", "v-vel"],
            timestep,
            datetime_str,
            use_cache,
        )

    reference_solution = zeta_solutions[0] if zeta_solutions else velocity_solutions[0]
    time_label = _format_time_label(reference_solution["time"])
    distance = reference_solution["distance"]

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

    if zeta_solutions:
        depth_positive_up = -np.array(zeta_solutions[0]["depth"])
        zeta_series = [
            {
                "values": solution["zeta"],
                "label": f"{solution_labels[i]} (zeta)",
            }
            for i, solution in enumerate(zeta_solutions)
        ]
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        plot_solution_along_transect(
            fig,
            ax,
            distance,
            zeta_series,
            ylabel=f"Elevation ({elev_unit})",
            title=f"Water Level and Bathymetry at {time_label}\n{transect_suffix}",
            plot_in_foot=args.plot_in_foot,
            ymin=args.ymin,
            ymax=args.ymax,
            solution_colors=solution_colors,
            secondary_values=depth_positive_up,
            secondary_label="Bathymetry (-depth)",
            secondary_color="saddlebrown",
        )
        _save_figure(fig, f"{args.output_prefix}_zeta.png")

    if velocity_solutions:
        velocity_plots = [
            ("u-vel", "u-vel", f"U-Velocity ({vel_unit})", "U-Velocity"),
            ("v-vel", "v-vel", f"V-Velocity ({vel_unit})", "V-Velocity"),
            (
                "v-tangential",
                "v-tangential",
                f"Tangential Velocity ({vel_unit})",
                "Tangential Velocity",
            ),
        ]

        for suffix, var_key, ylabel, plot_name in velocity_plots:
            series = []
            for i, solution in enumerate(velocity_solutions):
                if var_key == "v-tangential":
                    values = solution["u-vel"] * tx + solution["v-vel"] * ty
                else:
                    values = solution[var_key]

                series.append(
                    {
                        "values": values,
                        "label": f"{solution_labels[i]} ({plot_name})",
                    }
                )

            fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
            plot_solution_along_transect(
                fig,
                ax,
                distance,
                series,
                ylabel=ylabel,
                title=f"{plot_name} at {time_label}\n{transect_suffix}",
                plot_in_foot=args.plot_in_foot,
                ymin=args.ymin,
                ymax=args.ymax,
                solution_colors=solution_colors,
            )
            _save_figure(fig, f"{args.output_prefix}_{suffix}.png")

    if zeta_solutions and velocity_solutions:
        flux_series = []
        for i, (zeta_solution, velocity_solution) in enumerate(
            zip(zeta_solutions, velocity_solutions)
        ):
            v_tangential = velocity_solution["u-vel"] * tx + velocity_solution["v-vel"] * ty
            total_depth = np.array(zeta_solution["zeta"]) + np.array(zeta_solution["depth"])
            flux_series.append(
                {
                    "values": total_depth * v_tangential,
                    "label": f"{solution_labels[i]} (flux)",
                }
            )

        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        flux_scale = (1.0 / 0.3048) ** 2 if args.plot_in_foot else 1.0
        plot_solution_along_transect(
            fig,
            ax,
            distance,
            flux_series,
            ylabel=f"Flux ({flux_unit})",
            title=f"Flux at {time_label}\n{transect_suffix}",
            unit_scale=flux_scale,
            ymin=args.ymin,
            ymax=args.ymax,
            solution_colors=solution_colors,
        )
        _save_figure(fig, f"{args.output_prefix}_flux.png")


if __name__ == "__main__":
    main()
