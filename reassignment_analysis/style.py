from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# Plotting order and colours
WR_ORDER = ['AM','COL','WCT','FH','TH','EH','CH','WH']  # Order of plotting WR
WR_COLORS = {  # Colors for WR
    'AM'  : "#b35c44",  # Terracotta
    'COL' : "#d39c55",  # Warm ochre
    'WCT' : "#a3a847",  # Olive green
    'FH'  : "#669966",  # Soft forest green
    'TH'  : "#5b8ea3",  # Dusty blue
    'EH'  : "#7b6ea7",  # Muted violet
    'CH'  : "#a87c7c",  # Earthy rose
    'WH'  : "#84786f",  # Weathered taupe
}

# General colors and sizes
BACKGROUND_COLOR = "#e8e6e1"
TEXT_COLOR = "#333333"
GRID_COLOR = "#444444"
LAND_COLOR = "lightgrey"
MATRIX_FIGSIZE = np.array([11.7/0.85, 8.3])
POPUP_FIGSIZE = MATRIX_FIGSIZE * 0.75
NVEC = 2  # Coarsening for quiver plots
DEFAULT_VARIABLES = ["msl","tcwv","uv_850hPa_TW","jet","pv_315K","pv_330K"] # Variables included in default plot

# --- Create linear colormaps for WRs ---
def make_linear_cmap(color, base=BACKGROUND_COLOR, name="custom_cmap"):
    return LinearSegmentedColormap.from_list(name, [base, color])

WR_CMAPS = {
    name: make_linear_cmap(color, BACKGROUND_COLOR, f"WR{name}_cmap")
    for name, color in WR_COLORS.items()
}

# --- Plot options for composite variables ---
VARIABLE_PLOTMODES = {
    "jet"         : ["contourf"],
    "tcwv"        : ["contourf", "contour"],
    "msl"         : ["contour"],
    "pv_315K"     : ["contour"],
    "pv_330K"     : ["contour"],
    "z_500hPa"    : ["contour"],
    "uv_850hPa_TW": ["quiver"],
    "uv_850hPa"   : ["quiver"],
    "uv_500hPa"   : ["quiver"]
}

VARIABLE_CMAPS = {
    "jet"         : "Blues",
    "tcwv"        : "Greens"
}

VARIABLE_ALPHAS = {
    "jet"         : 0.75,
    "tcwv"        : 0.4,
    "gridlines"   : 0.5
}

VARIABLE_CFLEVELS = {
    "jet"         : np.arange(20, 52.5, 2.5),
    "tcwv"        : np.arange(48, 1000, 500)
}

VARIABLE_COLORS = {
    "msl"         : "black",
    "tcwv"        : "limegreen",
    "z_500hPa"    : "C0",
    "pv_330K"     : "magenta",
    "pv_315K"     : "magenta",
    "uv_850hPa_TW": "darkgreen",
    "uv_850hPa"   : "black",
    "uv_500hPa"   : "black",
    "gridlines"   : "lightgrey"
}

VARIABLE_LINESTYLES = {
    "msl"         : "-",
    "tcwv"        : "-",
    "z_500hPa"    : "-",
    "pv_330K"     : "--",
    "pv_315K"     : "-",
    "gridlines"   : "--"
}

VARIABLE_LINEWIDTHS = {
    "msl"         : 0.75,
    "tcwv"        : 0.8,
    "z_500hPa"    : 0.75,
    "pv_330K"     : 1.0,
    "pv_315K"     : 1.25,
    "coastline"   : 0.4,
    "gridlines"   : 0.4
}

VARIABLE_CLEVELS = {
    "msl"         : np.arange(900, 1064, 4),
    "tcwv"        : [45],
    "z_500hPa"    : np.arange(-5.1, 5.1, 0.2),
    "pv_330K"     : [-2],
    "pv_315K"     : [-2]
}
# Add clickable tickboxes for plot options
LEGEND_LABELS = {
  'msl': r'$\overline{msl}$',
  'tcwv': r'$\overline{tcwv}$',
  'jet': 'jet',
  'pv_315K': r'$\overline{pv}$ 315K',
  'pv_330K': r'$\overline{pv}$ 330K',
  'z_500hPa': r"$\overline{z'}$ 500hPa",
  'uv_850hPa_TW': r'$\overline{u},\overline{v}$ 850hPa TW',
  'uv_850hPa': r'$\overline{u},\overline{v}$ 850hPa',
  'uv_500hPa': r'$\overline{u},\overline{v}$ 500hPa'
}

QUIVER_SETTINGS = {
  "scale" : 3,
  "width" : 0.003,
  "minshaft" = 2
}
