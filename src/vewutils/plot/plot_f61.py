"""Plot fort.61.nc water levels for selected stations on one axes."""

from __future__ import annotations

import sys
from datetime import datetime, timedelta
from pathlib import Path


def parse_station_ranges(tokens: list[str]) -> list[tuple[int, int]]:
    """Parse station tokens into inclusive ``(start, end)`` appearance orders.

    Orders are 1-based. ``1-4,8`` and ``1 5 12`` are both accepted. Ranges are
    not expanded here.
    """
    ranges: list[tuple[int, int]] = []
    for token in tokens:
        for part in token.split(','):
            part = part.strip()
            if not part:
                continue
            bounds = part.split('-')
            if len(bounds) == 1:
                start = end = _parse_positive_order(bounds[0], part)
            elif len(bounds) == 2 and bounds[0] and bounds[1]:
                start = _parse_positive_order(bounds[0], part)
                end = _parse_positive_order(bounds[1], part)
                if end < start:
                    raise ValueError(
                        f'station range {part!r} is descending; '
                        f'write it as {end}-{start}'
                    )
            else:
                raise ValueError(
                    f'invalid station order {part!r}; '
                    f'use a positive integer or an inclusive range such as 1-4'
                )
            ranges.append((start, end))
    if not ranges:
        raise ValueError(
            'no station orders given; pass 1-based appearance orders '
            'such as 1 5 12 or 1-4,8'
        )
    return ranges


def expand_station_orders(
        ranges: list[tuple[int, int]],
        nstation: int) -> list[int]:
    """Expand inclusive ranges to appearance orders within ``1..nstation``."""
    if nstation < 1:
        raise ValueError('fort.61.nc station dimension is empty')
    orders: list[int] = []
    seen: set[int] = set()
    for start, end in ranges:
        if start > nstation or end > nstation:
            outside = end if end > nstation else start
            raise ValueError(
                f'station order {outside} is outside 1..{nstation}'
            )
        for order in range(start, end + 1):
            if order in seen:
                raise ValueError(f'duplicate station order {order}')
            seen.add(order)
            orders.append(order)
    return orders


def _parse_positive_order(text: str, token: str) -> int:
    try:
        order = int(text)
    except ValueError as exc:
        raise ValueError(
            f'invalid station order {token!r}; '
            f'use a positive integer or an inclusive range such as 1-4'
        ) from exc
    if order < 1:
        raise ValueError(
            f'station order {token!r} must be a positive 1-based appearance order'
        )
    return order


def load_f61_stations(path: str | Path, tokens: list[str]):
    """Read water levels for the stations named by ``tokens``.

    Returns
    -------
    orders : list of int
        1-based appearance orders, in the requested order.
    times : pandas.DatetimeIndex
        Model times, timezone-naive.
    water_levels : ndarray
        Shape ``(len(orders), ntime)``, with fill values replaced by NaN.
    """
    import numpy as np
    import xarray as xr

    ranges = parse_station_ranges(tokens)
    path = Path(path)
    with xr.open_dataset(path) as ds:
        if 'zeta' not in ds:
            raise ValueError(f'{path}: no zeta variable')
        zeta = ds['zeta']
        if 'station' not in zeta.dims:
            raise ValueError(f'{path}: zeta has no station dimension')
        nstation = int(zeta.sizes['station'])
        orders = expand_station_orders(ranges, nstation)
        selected = zeta.isel(station=[order - 1 for order in orders])
        water_levels = np.asarray(selected.values, dtype=float)
        station_axis = selected.dims.index('station')
        water_levels = np.moveaxis(water_levels, station_axis, 0)
        fill_value = _zeta_fill_value(zeta)
        if fill_value is not None and np.isfinite(fill_value):
            water_levels = np.where(
                water_levels == fill_value, np.nan, water_levels
            )
        times = _time_index_from_dataset(ds, path)
    if water_levels.shape[1] != len(times):
        raise ValueError(
            f'{path}: zeta time length {water_levels.shape[1]} '
            f'does not match time variable length {len(times)}'
        )
    return orders, times, water_levels


def _zeta_fill_value(zeta):
    raw = zeta.attrs.get('_FillValue', zeta.encoding.get('_FillValue'))
    if raw is None:
        raw = zeta.attrs.get('missing_value', -99999.0)
    import numpy as np

    try:
        fill_value = float(np.asarray(raw).reshape(-1)[0])
    except (TypeError, ValueError):
        return None
    return fill_value


def _time_index_from_dataset(ds, path: Path):
    """Return a timezone-naive DatetimeIndex for ``ds['time']``."""
    import numpy as np
    import pandas as pd

    if 'time' not in ds:
        raise ValueError(f'{path}: no time variable')
    tvar = ds['time']
    values = tvar.values
    if len(np.asarray(values)) == 0:
        raise ValueError(f'{path}: time dimension is empty')
    if np.issubdtype(np.asarray(values).dtype, np.number):
        units = str(tvar.attrs.get('units', ''))
        if 'since' in units:
            origin = units.split('since', 1)[1].strip()
            times = pd.to_datetime(values, unit='s', origin=origin)
        else:
            base_date = tvar.attrs.get('base_date')
            if base_date is None:
                raise ValueError(
                    f'{path}: numeric time values without decodable units/base_date'
                )
            origin = pd.Timestamp(str(base_date))
            times = origin + pd.to_timedelta(values, unit='s')
    else:
        times = pd.to_datetime(np.asarray(values).astype('datetime64[ms]'))
    times = pd.DatetimeIndex(times)
    if times.tz is not None:
        times = times.tz_convert('UTC').tz_localize(None)
    return times


def resolve_xlim(times, date_start_str: str | None, date_end_str: str | None):
    """Return plot limits, or ``(None, None)`` when both dates are omitted.

    A provided ``date-end`` is an inclusive calendar day, so the right limit is
    the following midnight.
    """
    if date_start_str is None and date_end_str is None:
        return None, None

    import pandas as pd

    try:
        start_day = (
            datetime.strptime(date_start_str, '%Y-%m-%d')
            if date_start_str is not None else None
        )
        end_day = (
            datetime.strptime(date_end_str, '%Y-%m-%d')
            if date_end_str is not None else None
        )
    except ValueError as exc:
        raise ValueError('dates must be YYYY-MM-DD') from exc

    if start_day is not None and end_day is not None and start_day > end_day:
        raise ValueError(
            f'date-start {date_start_str} is after date-end {date_end_str}'
        )

    t_min = pd.Timestamp(times.min()).to_pydatetime().replace(tzinfo=None)
    t_max = pd.Timestamp(times.max()).to_pydatetime().replace(tzinfo=None)
    left = t_min if start_day is None else start_day
    # date-end is an inclusive calendar day, so the right edge is the next midnight
    right = t_max if end_day is None else end_day + timedelta(days=1)
    if left >= right:
        raise ValueError(
            f'plot window {left} .. {right} does not overlap the file '
            f'({t_min} .. {t_max})'
        )
    return left, right


def plot_f61(
        ax, times, water_levels, labels,
        date_start=None, date_end=None,
        colors=None, linestyles=None):
    """Overlay one water-level series per station on ``ax``.

    ``colors`` and ``linestyles`` are one entry per series. ``None`` keeps the
    matplotlib default for that property (color cycle, solid lines).
    """
    from matplotlib.dates import DateFormatter

    nseries = len(water_levels)
    if len(labels) != nseries:
        raise ValueError(f'expected {nseries} labels, got {len(labels)}')
    if colors is not None and len(colors) != nseries:
        raise ValueError(f'expected {nseries} colors, got {len(colors)}')
    if linestyles is not None and len(linestyles) != nseries:
        raise ValueError(
            f'expected {nseries} linestyles, got {len(linestyles)}'
        )
    for i, (series, label) in enumerate(zip(water_levels, labels)):
        kwargs = {'label': label}
        if colors is not None:
            kwargs['color'] = colors[i]
        if linestyles is not None:
            kwargs['linestyle'] = linestyles[i]
        ax.plot(times, series, **kwargs)
    ax.set_ylabel('Water Level (m)')
    ax.grid(True)
    ax.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M'))
    ax.legend(loc='best')
    if date_start is not None or date_end is not None:
        ax.set_xlim(date_start, date_end)


def legend_labels(orders: list[int], labels: list[str] | None) -> list[str]:
    """Return legend text: appearance orders, or ``labels`` when provided."""
    if labels is None:
        return [str(order) for order in orders]
    return _match_series_values(labels, len(orders), '--labels')


def series_colors(nseries: int, colors: list[str] | None) -> list[str] | None:
    """Return one color per series, or ``None`` to use the default color cycle."""
    if colors is None:
        return None
    parsed = _match_series_values(_split_list_tokens(colors), nseries, '--colors')
    _validate_colors(parsed)
    return parsed


def series_linestyles(
        nseries: int, linestyles: list[str] | None) -> list[str] | None:
    """Return one linestyle per series, or ``None`` to draw solid lines."""
    if linestyles is None:
        return None
    parsed = _match_series_values(
        _split_list_tokens(linestyles), nseries, '--linestyles'
    )
    _validate_linestyles(parsed)
    return parsed


def _split_list_tokens(tokens: list[str]) -> list[str]:
    """Split argparse tokens on commas. ``b,r,g`` and ``b r g`` both work."""
    items: list[str] = []
    for token in tokens:
        for part in token.split(','):
            part = part.strip()
            if part:
                items.append(part)
    return items


def _match_series_values(values: list[str], nseries: int, option: str) -> list[str]:
    if len(values) != nseries:
        raise ValueError(
            f'{option} has {len(values)} entries for {nseries} stations; '
            f'pass one value per station in expanded order'
        )
    return list(values)


def _validate_colors(colors: list[str]) -> None:
    from matplotlib.colors import to_rgba

    for color in colors:
        try:
            to_rgba(color)
        except ValueError as exc:
            raise ValueError(f'invalid line color {color!r}') from exc


def _validate_linestyles(linestyles: list[str]) -> None:
    from matplotlib.lines import lineStyles, ls_mapper

    accepted = set(lineStyles) | set(ls_mapper) | set(ls_mapper.values())
    valid = 'solid, dashed, dotted, dashdot, -, --, -., :'
    for style in linestyles:
        if style not in accepted:
            raise ValueError(
                f'invalid line style {style!r}. Valid styles: {valid}'
            )


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        add_help=False,
        description=(
            'Plot modeled water-level hydrographs from a fort.61.nc file on '
            'one axes. Stations are selected by 1-based appearance order along '
            'the station dimension. No observations are plotted.'
        ),
    )
    parser.add_argument(
        '--f61',
        required=True,
        help='Path to a fort.61.nc file',
    )
    parser.add_argument(
        '--stations',
        required=True,
        nargs='+',
        metavar='ORDER',
        help=(
            '1-based station appearance orders. Commas and inclusive ranges '
            'are expanded (for example: 1 5 12 or 1-4,8).'
        ),
    )
    parser.add_argument(
        '--labels',
        nargs='+',
        default=None,
        help=(
            'Legend labels in the same order as the expanded --stations list. '
            'Default: the appearance-order integers.'
        ),
    )
    parser.add_argument(
        '--colors',
        nargs='+',
        default=None,
        help=(
            'Line colors in the same order as the expanded --stations list '
            '(for example: b r g, or --colors=b,r,g). '
            'Default: matplotlib color cycle.'
        ),
    )
    parser.add_argument(
        '--linestyles',
        nargs='+',
        default=None,
        help=(
            'Line styles in the same order as the expanded --stations list. '
            'Names: solid, dashed, dotted, dashdot. Short codes that begin '
            'with a dash must be one comma-separated value, for example '
            '--linestyles=-,--,:. Default: solid.'
        ),
    )
    parser.add_argument(
        '--date-start',
        default=None,
        help=(
            'Inclusive plot start date (YYYY-MM-DD). '
            'Default: first time in the file.'
        ),
    )
    parser.add_argument(
        '--date-end',
        default=None,
        help=(
            'Inclusive plot end date (YYYY-MM-DD). '
            'Default: last time in the file.'
        ),
    )
    parser.add_argument(
        '--outputfile',
        required=True,
        help='Output figure file name',
    )
    parser.add_argument(
        '--fig-width',
        type=float,
        default=10.0,
        help='Figure width in inches (default: 10.0)',
    )
    parser.add_argument(
        '--fig-height',
        type=float,
        default=6.0,
        help='Figure height in inches (default: 6.0)',
    )
    return parser


def main(args=None):
    import matplotlib.pyplot as plt

    if args is None:
        args = get_parser().parse_args()

    try:
        orders, times, water_levels = load_f61_stations(args.f61, args.stations)
        labels = legend_labels(orders, args.labels)
        colors = series_colors(len(orders), args.colors)
        linestyles = series_linestyles(len(orders), args.linestyles)
        date_start, date_end = resolve_xlim(
            times, args.date_start, args.date_end
        )
    except (OSError, ValueError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        return 1

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    plot_f61(
        ax, times, water_levels, labels, date_start, date_end,
        colors=colors, linestyles=linestyles,
    )
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(args.outputfile)
    plt.close(fig)
    print(f'Wrote {args.outputfile}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
