"""Named unit conversions applicable to plotted field values.

Used by plot_max_ele_2d.py and plot_solution_2d.py's `conversion` option, so
a report can plot a field in different units than however the source NetCDF
stores it (e.g. ADCIRC's pressure_min in equivalent meters of water, plotted
in hPa instead).
"""

from __future__ import annotations

# Each entry is (scale, offset), applied as `value * scale + offset`.
FIELD_CONVERSIONS: dict[str, tuple[float, float]] = {
    # ADCIRC stores atmospheric pressure as an equivalent water-column height
    # (see src/wind.F, src/constants.F90): P_mH2O = 100 * P_hPa / (rho0 * g),
    # with rho0 = 1000 kg/m^3 and g = 9.80665 m/s^2. Inverting:
    # P_hPa = P_mH2O * rho0 * g / 100 = P_mH2O * 98.0665.
    'mwater_to_hpa': (98.0665, 0.0),
    'm_to_ft': (3.28084, 0.0),
}


def apply_conversion(var_data, conversion: str | None):
    """Apply a named conversion (see FIELD_CONVERSIONS) to var_data.

    Returns var_data unchanged if conversion is None.
    """
    if conversion is None:
        return var_data
    try:
        scale, offset = FIELD_CONVERSIONS[conversion]
    except KeyError:
        valid = ', '.join(sorted(FIELD_CONVERSIONS))
        raise ValueError(
            f'Unknown conversion {conversion!r}. Available: {valid}'
        ) from None
    return var_data * scale + offset
