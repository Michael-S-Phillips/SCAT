import tkinter as tk
from tkinter import filedialog, messagebox, ttk, Button, Entry, Toplevel, Listbox, END
import pickle
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import spectral
import numpy as np
import rasterio as rio
from rasterio.warp import transform as rio_transform
import cv2
import geopandas as gpd
from shapely.geometry import Polygon as sgp
from scipy import interpolate
import pandas as pd
import glob
import re
import warnings
import matplotlib.colors as mcolors
from joblib import Parallel, delayed
import os

# Import the new plotting system
from plot_manager import get_plot_manager
from spectral_plot_window import SpectralPlotWindow, SpectrumData
from drag_drop_manager import get_drag_drop_manager

# try:
from hypyrameter.paramCalculator import cubeParamCalculator
from hypyrameter.interpNans import interpNaNs as interpNaNs
from hypyrameter.utils import (
    getSpecFiles,
    getWavelengthFromUSGS,
    getReflectanceFromUSGS,
)

# except:
#     print("Unable to load HyPyRameter, visit https://github.com/Michael-S-Phillips/HyPyRameter for more information.")


def _nanmean_over_mask(data, mask, fill_values=(), valid_range=None):
    """Per-band NaN-aware mean over pixels selected by a 2D mask.

    Replacement for spectral.calc_stats(...).mean, whose underlying
    np.average propagates NaN — a single NaN pixel poisons the whole
    band's mean. This uses np.nanmean so an ROI containing some bad
    pixels still produces a usable spectrum.

    Pre-filter (per pixel, per band, applied before averaging):
      - NaN is always treated as missing (np.nanmean).
      - Any value equal to one of `fill_values` (e.g. CRISM's 65535
        data-ignore sentinel from the ENVI header) is treated as NaN.
      - Values outside `valid_range = (lo, hi)` (either bound may be
        None) are treated as NaN.

    Filtering per-pixel — rather than computing the mean and then
    NaN-ing out-of-range *band means* — is what lets an ROI that
    straddles a no-data area still produce a real spectrum: the
    sentinel pixels are excluded one band at a time, so any band
    that has at least one valid pixel under the ROI gets a real
    mean.

    Returns a length-`bands` array. Bands where every selected pixel
    is NaN/sentinel/out-of-range return NaN (warning suppressed).
    """
    mask_bool = np.asarray(mask).astype(bool)
    # SPy's ImageArray (from .load()) overrides __getitem__ in a way that
    # mishandles 2D boolean fancy indexing — it iterates the mask rows and
    # raises "too many indices for array". View as a plain ndarray first.
    arr = np.asarray(data)
    if arr.ndim == 3:
        # SPy's load() returns (rows, cols, bands).
        pixels = arr[mask_bool]
        n_bands = arr.shape[-1]
    elif arr.ndim == 2:
        pixels = arr[mask_bool]
        n_bands = 1
    else:
        raise ValueError(f"Unsupported data ndim={arr.ndim}")

    if pixels.size == 0:
        return np.full(n_bands, np.nan, dtype=np.float64)

    pixels = pixels.astype(np.float64, copy=True)
    for fv in fill_values:
        # NaN sentinels need np.isnan; finite sentinels use ==.
        if fv is None:
            continue
        if isinstance(fv, float) and np.isnan(fv):
            continue  # already handled by nanmean
        pixels[pixels == fv] = np.nan
    if valid_range is not None:
        lo, hi = valid_range
        if lo is not None:
            pixels[pixels < lo] = np.nan
        if hi is not None:
            pixels[pixels > hi] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(pixels, axis=0)


def _data_ignore_value(spy_image):
    """Pull `data ignore value` from an ENVI header if present, as float.

    Returns None if the image isn't an SPy SpyFile-like object or the
    header doesn't set one. CRISM tiles set this to 65535.
    """
    md = getattr(spy_image, 'metadata', None)
    if not md:
        return None
    raw = md.get('data ignore value')
    if raw is None:
        return None
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


# Patterns describing how to derive the geometry/backplane file path
# from a CRISM source filename. Each entry is (regex, replacement); the
# first match wins. Order matters — `mrral` (mosaic tile albedo) must
# come before generic per-strip patterns because the per-strip regex
# would not match `mrral` anyway, but listing tile patterns first
# makes the intent explicit.
_GEOMETRY_NAME_PATTERNS = (
    # Map-projected mosaic tile reflectance products → DDR companion.
    # e.g. t0886_mrral_05s058_0327_4 → t0886_mrrde_05s058_0327_4
    (re.compile(r'_mrr(al|if|ir|sr|su|ra)_', re.IGNORECASE), '_mrrde_'),
    # Per-strip MTRDR products → input-geometry companion.
    # e.g. frt00003e12_07_if165j_mtr3 → frt00003e12_07_in165j_mtr3
    (re.compile(r'_(if|sr|su)(\d{3}[a-z])', re.IGNORECASE), r'_in\2'),
)


def _resolve_geometry_file(source_filename):
    """Map a CRISM source `.hdr`/`.img` path to its geometry companion.

    Returns the absolute path of the geometry `.img` file if a known
    naming convention matches AND the file exists on disk; otherwise
    None. The choice of companion depends on product type:

      - `*_mrr??_*` (map-projected mosaic tiles) → `*_mrrde_*`
      - `*_if|sr|suNNNx*`  (per-strip MTRDR)    → `*_inNNNx*`

    Returning None is the unambiguous "no geometry available" signal —
    callers should disable column-lock for this image rather than
    silently falling back to reading the source as a fake backplane.
    """
    directory = os.path.dirname(source_filename)
    filename = os.path.basename(source_filename)
    stem, _ = os.path.splitext(filename)

    for pattern, replacement in _GEOMETRY_NAME_PATTERNS:
        new_stem, n_subs = pattern.subn(replacement, stem)
        if n_subs == 0 or new_stem == stem:
            continue
        candidate = os.path.join(directory, new_stem + '.img')
        if os.path.exists(candidate):
            return candidate
    return None


# Band descriptions we look up in CRISM geometry files. The same logical
# layer is named differently across product generations (per-strip _in
# files vs tile _mrrde files), so we accept any of several aliases and
# match case-insensitively against rasterio's description strings.
#
# Note on strip identification for tiles: `target_id` is the actual
# unique strip identifier (one CRISM observation ID per source strip),
# *not* `segment_id`. Segment ID is a within-target counter that takes
# only ~4 values across an entire mosaic, so many distinct strips
# share the same Segment ID. Use Target ID for primary gating and
# Segment ID only as a secondary gate within multi-segment targets.
_BAND_ALIASES = {
    'ir_sample': (
        'ir (l-detector) sample',
        'ir sample',
        'l-detector sample',
    ),
    'target_id': (
        'target id',
    ),
    'segment_id': (
        'segment id (counter)',
        'segment id',
        'segment',
    ),
}


def _find_band_index(rio_image, kind):
    """Return 1-indexed rasterio band index for a logical layer.

    `kind` is a key into `_BAND_ALIASES`. Matches are case-insensitive
    and substring-tolerant (rasterio descriptions vary in punctuation
    and parenthetical detail across product versions). Returns None if
    no band description matches — caller must handle absence
    explicitly (e.g. tile DDR has Segment ID; per-strip _in files
    may not).
    """
    aliases = _BAND_ALIASES.get(kind, ())
    descs = rio_image.descriptions or ()
    for idx, desc in enumerate(descs, start=1):
        if not desc:
            continue
        d = desc.strip().lower()
        for alias in aliases:
            if alias in d:
                return idx
    return None


class SpectralCubeAnalysisTool:
    """
    TODO:
        - Adding more functionality to the polygon shape layer tool:
            > naming polygons
            > ability to edit individual polygons
            > plot spectrum from individual polygon
            > plot mean of multiple polygons together
            > ability to cancel a polygon in the middle of drawing
            > ratio polygons, duplicate polygons for ratioing, duplicate polygons in-column
        - Spectral plotting
            > add ability to save plots at figure quality
            > spectral smoothing routines (Rogers for CRISM, generic smoothing)
            > add ability to plot multiple spectra on the same plot for point and click spectra
            > add ability to plot library spectra
            > add ability to offset for clarity
        - Adding basic spectral analysis workflows
            > MNF
            > HyPyRameter
            > Vectorization (Sam)
            > Neural Net classifiers
        - Add default Browse product combos for parameter display
        - Add ability to save statistics from parameter file over ROIs
        - Add a tool for measuring distance.
        - Test with larger images.
        - Test with tripod images (will break because georef info won't be present)
        - Overall better error handling
        - adding interactive bad band list selection
            > show the image, be able to scroll through bands visually and also plot spectra,
              then use a span selector to select ranges that are bad bands.
        - abitility to drag and drop histogram lines
        - preset RGB combos for band parameter images
        - undo bad band select

    This class creates a GUI for analyzing hyperspectral data. It allows the user to load a hyperspectral and analyze an image
    """

    # ----------------------------------------------------------------
    # initial setup
    # ----------------------------------------------------------------
    def __init__(self, root):
        """
        Initializes the GUI and creates the main UI elements.
        """
        self.root = root
        self.root.title("Spectral Cube Analysis Tool")

        self.default_rgb_bands = [
            3,
            2,
            1,
        ]  # Fallback indices for RGB display of spectral data
        # Target wavelengths for CRISM spectral cubes (in nm)
        self.default_rgb_wavelengths = [2330.0, 1500.0, 774.0]  # R, G, B
        self.default_parameter_bands = [
            3,
            2,
            1,
        ]  # Fallback indices for RGB display of parameter data
        # Target band names for CRISM parameter images
        self.default_parameter_names = ["BD1300", "RPEAK1", "LCPINDEX2"]  # R, G, B
        self.default_stretch = "Percentile"  # Default stretch method for RGB display
        self.create_main_ui()
        self.create_menu()
        self.spectral_window = None
        self.ratio_spectral_window = None
        self.polygons_spectral_window = None
        self.polygons_ratio_spectral_window = None
        self.template_index_cache = None
        self.template_index = None
        self.denominator_index = None
        self.right_hist_window = None
        self.left_hist_window = None
        self.right_is_parameter = False
        self.left_is_parameter = False
        self.reset_right_hist = False
        self.reset_left_hist = False
        self.draw_polygons = False
        self.ignore_bad_bands_flag = False
        self.ignore_polygons_bad_bands_flag = False
        self.ignore_ratio_bad_bands_flag = False
        self.draw_poly_motion_flag = False
        self.lock_column = False
        # Geometry/backplane cache for column-lock. Keyed by source
        # filename. Value is a dict (when geometry was resolved and
        # loaded), or None (cached "no geometry available" sentinel
        # so we don't re-attempt resolution on every paste). Dict
        # entry: {'ir_sample': np.ndarray, 'segment_id': np.ndarray
        # | None, 'nodata': float | None, 'path': str}.
        self._geometry_cache = {}
        # Source files for which we've already shown the user the
        # "no geometry, column-lock disabled" warning — keep that
        # one-shot per cube so we don't spam dialogs.
        self._geometry_warning_shown = set()
        # Strip-footprint overlay state. The overlay outlines the
        # source CRISM strip (target_id + segment_id) on both canvases
        # while column-lock is on, so users can see the y-range and
        # x-range within which their next paste can actually lock.
        # Artist refs are kept so we can remove them on toggle/redraw.
        self._strip_footprint_artists = []
        # Cache of contour polygons keyed by (target_id, segment_id).
        # Recomputing the boundary of a 1.6M-pixel mask isn't cheap
        # enough to do on every redraw.
        self._strip_footprint_cache = {}
        self.points = []
        self.current_polygon = []
        self.denom_polygon = []
        self.polygon_table_item_id_row = []
        self.polygon_to_highlight = False
        self.all_polygons = []
        self.polygon_colors = []
        self.polygon_spectra = []
        self.polygons_library_reflectance = []
        self.polygon_params = []
        self.polygon_ratio_spectra = []
        self.polygons_ratio_library_reflectance = []
        self.spectrum = []
        self.denom_spectrum = []
        self.ratio_spectrum = []
        self.ratio_spectrum_list = []
        self.num_idx = None
        self.best_denom_id_list = []
        self.add_spectrum_flag = False
        self.bad_bands = None
        # Library spectra tracking for stretch/offset controls
        # Format: {name: {'line': line_obj, 'original_y': array, 'stretch': 1.0, 'offset': 0.0, 'anchor': value, 'pivot': value, 'wvl': array}}
        self.library_spectra_data = {}
        self.library_spectra_controls_frame = None
        self.collected_points = pd.DataFrame(
            columns=(
                "name",
                "color",
                "mineral_id1",
                "mineral_id2",
                "mineral_id3",
                "mineral_id4",
                "x",
                "y",
                "lon",
                "lat",
                "spectrum",
                "denom",
                "parameter_values",
            )
        )
        self.polygons_on_flag = False
        self.upload_polygons_flag = False

        self.usgs_spectra_path = "librarySpectra/"
        self.usgs_spectra_folders = glob.glob(self.usgs_spectra_path + "/*")
        self.mica_spectra_path = "mrocr_8001/mrocr_8001/data/"
        self.mica_spectra_names = glob.glob(self.mica_spectra_path + "*.tab")

        self.pc = None

        # Initialize the new plotting system
        self.plot_manager = get_plot_manager()
        self.drag_drop_manager = get_drag_drop_manager(self.root)

    def create_main_ui(self):
        """
        Creates the main UI elements for the application.
        """
        # Create a frame to hold UI elements with a fixed size
        self.main_ui_frame = tk.Frame(self.root)
        self.main_ui_frame.grid(
            row=0, column=1, sticky=tk.N + tk.S + tk.W + tk.E
        )  # Use grid and place it in column 1
        self.root.grid_rowconfigure(
            0, weight=1
        )  # Allow row 0 to expand vertically with window resize
        self.root.grid_columnconfigure(
            1, weight=1
        )  # Allow column 1 to expand horizontally with window resize

        # Create a label to display coordinates
        self.coordinates_label = tk.Label(
            self.main_ui_frame, text="Lat: 0.0000, Lon: 0.0000"
        )
        self.coordinates_label.grid(row=0, column=0, columnspan=2, sticky="nsew")
        # self.coordinates_label.grid_rowconfigure(row=0, weight=0)

        # frame for the buttons on the righthand side
        self.right_buttons_frame = tk.Frame(self.main_ui_frame)
        self.right_buttons_frame.grid(row=0, column=2, rowspan=2, sticky="new")
        self.main_ui_frame.grid_columnconfigure(2, weight=1)
        self.right_buttons_frame.grid_rowconfigure(0, weight=1)
        self.right_buttons_frame.grid_columnconfigure(0, weight=1)

        # Create a button to reset displays
        self.reset_display_button = tk.Button(
            self.right_buttons_frame,
            text="Reset Display Extent",
            command=self.reset_display_extent,
        )
        self.reset_display_button.grid(row=0, column=0, sticky="new")
        # self.reset_display_button.grid_rowconfigure(row=0, weight=0)

        # Add a button to draw polygons
        self.draw_polygons_button = tk.Button(
            self.right_buttons_frame,
            text="Draw Polygons",
            command=self.create_polygons_menu_window,
        )
        self.draw_polygons_button.grid(row=1, column=0, sticky="new")

        # Add a button to calculate spectral parameters
        # self.parameter_calculation_button = tk.Button(self.right_buttons_frame, text="Calculate Spectral Parameters", command=self.calculate_spectral_parameters)
        # self.parameter_calculation_button.grid(row=2, column=0, sticky = 'new')

        # Create a sub-frame for the left canvas buttons
        self.button_frame = tk.Frame(self.main_ui_frame)
        self.button_frame.grid(
            row=2, column=0, sticky=tk.N + tk.S + tk.W + tk.E
        )  # Use grid and place it in column 1
        self.main_ui_frame.grid_rowconfigure(
            0, weight=0
        )  # Allow row 0 to expand vertically with window resize
        self.main_ui_frame.grid_columnconfigure(
            1, weight=1
        )  # Allow column 1 to expand horizontally with window resize

        # Create band selection options for left canvas
        self.left_red_band_label = tk.Label(self.button_frame, text="Red Band:")
        self.left_red_band_label.grid(row=0, column=0, sticky=tk.W)

        self.left_red_band_var = tk.StringVar()
        self.left_red_band_menu = ttk.Combobox(
            self.button_frame, textvariable=self.left_red_band_var, state="readonly"
        )
        self.left_red_band_menu.grid(row=0, column=1, sticky=tk.W)

        self.left_green_band_label = tk.Label(self.button_frame, text="Green Band:")
        self.left_green_band_label.grid(row=1, column=0, sticky=tk.W)

        self.left_green_band_var = tk.StringVar()
        self.left_green_band_menu = ttk.Combobox(
            self.button_frame, textvariable=self.left_green_band_var, state="readonly"
        )
        self.left_green_band_menu.grid(row=1, column=1, sticky=tk.W)

        self.left_blue_band_label = tk.Label(self.button_frame, text="Blue Band:")
        self.left_blue_band_label.grid(row=2, column=0, sticky=tk.W)

        self.left_blue_band_var = tk.StringVar()
        self.left_blue_band_menu = ttk.Combobox(
            self.button_frame, textvariable=self.left_blue_band_var, state="readonly"
        )
        self.left_blue_band_menu.grid(row=2, column=1, sticky=tk.W)

        self.apply_button = tk.Button(
            self.button_frame, text="Apply", command=self.apply_new_left_bands
        )
        self.apply_button.grid(row=3, column=0, columnspan=2, sticky=tk.W)

        # Stretch type selector for left canvas
        self.left_stretch_type_label = tk.Label(self.button_frame, text="Stretch:")
        self.left_stretch_type_label.grid(row=4, column=0, sticky=tk.W)

        self.left_stretch_type_var = tk.StringVar(value=self.default_stretch)
        self.left_stretch_type_menu = ttk.Combobox(
            self.button_frame, textvariable=self.left_stretch_type_var,
            values=["Percentile", "Sigma Clip", "Hist. Equalization"],
            state="readonly", width=15
        )
        self.left_stretch_type_menu.grid(row=4, column=1, sticky=tk.W)

        # Create RGB stretch value options for left canvas
        self.left_red_stretch_label = tk.Label(self.button_frame, text="Red Stretch:")
        self.left_red_stretch_label.grid(row=5, column=0, sticky=tk.W)

        self.left_red_min_stretch_var = tk.DoubleVar(value=0.0)
        self.left_red_min_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.left_red_min_stretch_var
        )
        self.left_red_min_stretch_entry.grid(row=6, column=0, sticky=tk.W)

        self.left_red_max_stretch_var = tk.DoubleVar(value=1.0)
        self.left_red_max_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.left_red_max_stretch_var
        )
        self.left_red_max_stretch_entry.grid(row=6, column=1, sticky=tk.W)

        self.left_green_stretch_label = tk.Label(
            self.button_frame, text="Green Stretch:"
        )
        self.left_green_stretch_label.grid(row=7, column=0, sticky=tk.W)

        self.left_green_min_stretch_var = tk.DoubleVar(value=0.0)
        self.left_green_min_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.left_green_min_stretch_var
        )
        self.left_green_min_stretch_entry.grid(row=8, column=0, sticky=tk.W)

        self.left_green_max_stretch_var = tk.DoubleVar(value=1.0)
        self.left_green_max_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.left_green_max_stretch_var
        )
        self.left_green_max_stretch_entry.grid(row=8, column=1, sticky=tk.W)

        self.left_blue_stretch_label = tk.Label(self.button_frame, text="Blue Stretch:")
        self.left_blue_stretch_label.grid(row=9, column=0, sticky=tk.W)

        self.left_blue_min_stretch_var = tk.DoubleVar(value=0.0)
        self.left_blue_min_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.left_blue_min_stretch_var
        )
        self.left_blue_min_stretch_entry.grid(row=10, column=0, sticky=tk.W)

        self.left_blue_max_stretch_var = tk.DoubleVar(value=1.0)
        self.left_blue_max_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.left_blue_max_stretch_var
        )
        self.left_blue_max_stretch_entry.grid(row=10, column=1, sticky=tk.W)

        self.apply_stretch_button = tk.Button(
            self.button_frame, text="Apply Stretch", command=self.update_left_display
        )
        self.apply_stretch_button.grid(row=11, column=0, sticky=tk.W)

        self.auto_stretch_left_button = tk.Button(
            self.button_frame, text="Auto Stretch", command=self.auto_stretch_left
        )
        self.auto_stretch_left_button.grid(row=11, column=1, sticky=tk.W)

        # ----------------------------------------------------------------
        # Create a sub-frame for the right canvas buttons
        # ----------------------------------------------------------------
        self.button_frame = tk.Frame(self.main_ui_frame)
        self.button_frame.grid(
            row=2, column=1, sticky=tk.N + tk.S + tk.W + tk.E
        )  # Use grid and place it in column 1
        self.main_ui_frame.grid_rowconfigure(
            0, weight=1
        )  # Allow row 0 to expand vertically with window resize
        self.main_ui_frame.grid_columnconfigure(
            1, weight=1
        )  # Allow column 1 to expand horizontally with window resize

        # Create band selection options for the right canvas
        self.right_red_band_label = tk.Label(self.button_frame, text="Red Band:")
        self.right_red_band_label.grid(row=0, column=0, sticky=tk.W)

        self.right_red_band_var = tk.StringVar()
        self.right_red_band_menu = ttk.Combobox(
            self.button_frame, textvariable=self.right_red_band_var, state="readonly"
        )
        self.right_red_band_menu.grid(row=0, column=1, sticky=tk.W)

        self.right_green_band_label = tk.Label(self.button_frame, text="Green Band:")
        self.right_green_band_label.grid(row=1, column=0, sticky=tk.W)

        self.right_green_band_var = tk.StringVar()
        self.right_green_band_menu = ttk.Combobox(
            self.button_frame, textvariable=self.right_green_band_var, state="readonly"
        )
        self.right_green_band_menu.grid(row=1, column=1, sticky=tk.W)

        self.right_blue_band_label = tk.Label(self.button_frame, text="Blue Band:")
        self.right_blue_band_label.grid(row=2, column=0, sticky=tk.W)

        self.right_blue_band_var = tk.StringVar()
        self.right_blue_band_menu = ttk.Combobox(
            self.button_frame, textvariable=self.right_blue_band_var, state="readonly"
        )
        self.right_blue_band_menu.grid(row=2, column=1, sticky=tk.W)

        self.apply_button = tk.Button(
            self.button_frame, text="Apply", command=self.apply_new_right_bands
        )
        self.apply_button.grid(row=3, column=0, columnspan=2, sticky=tk.W)

        # Stretch type selector for right canvas
        self.right_stretch_type_label = tk.Label(self.button_frame, text="Stretch:")
        self.right_stretch_type_label.grid(row=4, column=0, sticky=tk.W)

        self.right_stretch_type_var = tk.StringVar(value=self.default_stretch)
        self.right_stretch_type_menu = ttk.Combobox(
            self.button_frame, textvariable=self.right_stretch_type_var,
            values=["Percentile", "Sigma Clip", "Hist. Equalization"],
            state="readonly", width=15
        )
        self.right_stretch_type_menu.grid(row=4, column=1, sticky=tk.W)

        # Create RGB stretch value options for the right canvas
        self.right_red_stretch_label = tk.Label(self.button_frame, text="Red Stretch:")
        self.right_red_stretch_label.grid(row=5, column=0, sticky=tk.W)

        self.right_red_min_stretch_var = tk.DoubleVar(value=0.0)
        self.right_red_min_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.right_red_min_stretch_var
        )
        self.right_red_min_stretch_entry.grid(row=6, column=0, sticky=tk.W)

        self.right_red_max_stretch_var = tk.DoubleVar(value=1.0)
        self.right_red_max_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.right_red_max_stretch_var
        )
        self.right_red_max_stretch_entry.grid(row=6, column=1, sticky=tk.W)

        self.right_green_stretch_label = tk.Label(
            self.button_frame, text="Green Stretch:"
        )
        self.right_green_stretch_label.grid(row=7, column=0, sticky=tk.W)

        self.right_green_min_stretch_var = tk.DoubleVar(value=0.0)
        self.right_green_min_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.right_green_min_stretch_var
        )
        self.right_green_min_stretch_entry.grid(row=8, column=0, sticky=tk.W)

        self.right_green_max_stretch_var = tk.DoubleVar(value=1.0)
        self.right_green_max_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.right_green_max_stretch_var
        )
        self.right_green_max_stretch_entry.grid(row=8, column=1, sticky=tk.W)

        self.right_blue_stretch_label = tk.Label(
            self.button_frame, text="Blue Stretch:"
        )
        self.right_blue_stretch_label.grid(row=9, column=0, sticky=tk.W)

        self.right_blue_min_stretch_var = tk.DoubleVar(value=0.0)
        self.right_blue_min_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.right_blue_min_stretch_var
        )
        self.right_blue_min_stretch_entry.grid(row=10, column=0, sticky=tk.W)

        self.right_blue_max_stretch_var = tk.DoubleVar(value=1.0)
        self.right_blue_max_stretch_entry = tk.Entry(
            self.button_frame, textvariable=self.right_blue_max_stretch_var
        )
        self.right_blue_max_stretch_entry.grid(row=10, column=1, sticky=tk.W)

        self.apply_stretch_button = tk.Button(
            self.button_frame, text="Apply Stretch", command=self.update_right_display
        )
        self.apply_stretch_button.grid(row=11, column=0, sticky=tk.W)

        self.auto_stretch_right_button = tk.Button(
            self.button_frame, text="Auto Stretch", command=self.auto_stretch_right
        )
        self.auto_stretch_right_button.grid(row=11, column=1, sticky=tk.W)

        # Configure rows and columns to expand with resizing
        self.button_frame.grid_columnconfigure(
            0, weight=1
        )  # Allow column 0 to expand horizontally with window resize
        self.button_frame.grid_columnconfigure(
            1, weight=1
        )  # Allow column 1 to expand horizontally with window resize

    def create_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(
            label="Load Hyperspectral Cube", command=self.load_left_data
        )
        file_menu.add_command(
            label="Load Band Parameter Image", command=self.load_right_data
        )
        file_menu.add_separator()
        file_menu.add_command(label="Load Session", command=self.load_state)
        file_menu.add_command(label="Save Session", command=self.save_state)
        file_menu.add_command(label="Exit", command=self.root.quit)

        histogram_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Plot Hitograms", menu=histogram_menu)
        histogram_menu.add_command(
            label="Plot Left Frame Histogram", command=self.plot_left_histograms
        )
        histogram_menu.add_command(
            label="Plot Right Frame Histogram", command=self.plot_right_histograms
        )

        # Plot windows menu
        plot_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Plot Windows", menu=plot_menu)
        plot_menu.add_command(
            label="New Spectral Plot Window", command=self.create_new_spectral_plot_window
        )

        processing_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Processing", menu=processing_menu)
        processing_menu.add_command(
            label="Calculate Spectral Parameters",
            command=self.calculate_spectral_parameters,
        )
        processing_menu.add_command(
            label="Select Bad Bands", command=self.select_bad_bands
        )

    def create_left_canvas(self):
        self.left_figure, self.left_ax = plt.subplots(figsize=(6, 6))
        self.left_ax.set_title(f"Hyperspectral Cube: {self.img_base_name}", fontsize=8)
        self.left_ax.tick_params(axis="both", which="major", labelsize=8)
        # Set the margins and spacing to zero
        self.left_figure.tight_layout()
        self.left_figure.subplots_adjust(wspace=0, hspace=0)

        self.left_frame = tk.Frame(self.main_ui_frame)
        self.left_frame.grid(row=1, column=0, sticky="nsew")

        self.left_canvas = FigureCanvasTkAgg(self.left_figure, master=self.left_frame)
        self.left_canvas.get_tk_widget().grid(
            row=0, column=0, sticky="nsew"
        )  # Use grid instead of pack

        self.left_nav_toolbar_frame = tk.Frame(self.left_frame)
        self.left_nav_toolbar = NavigationToolbar2Tk(
            self.left_canvas, self.left_nav_toolbar_frame
        )
        self.left_nav_toolbar.update()
        self.left_nav_toolbar_frame.grid(row=1, column=0, sticky="nsew")

        # Configure columns to expand with resizing
        self.left_frame.grid_rowconfigure(
            0, weight=1
        )  # Adjust the column number if needed
        self.left_frame.grid_columnconfigure(0, weight=1)
        self.left_canvas.get_tk_widget().grid_columnconfigure(0, weight=1)
        self.left_canvas.get_tk_widget().grid_rowconfigure(0, weight=1)

        if self.left_data.nrows / self.left_data.ncols > 2:
            # Create variables to store the current top-left pixel coordinates
            self.top_left_row = 0
            self.top_left_col = 0

            # Define the fixed size of the displayed portion
            self.left_display_rows = int(3 * self.left_data.ncols)
            self.left_ax.set_ylim(int(3 * self.left_data.ncols), 0)
            self.left_display_cols = self.left_data.ncols
            self.left_ax.set_xlim(0, self.left_data.ncols)
        else:
            # Create variables to store the current top-left pixel coordinates
            self.top_left_row = 0
            self.top_left_col = 0

            # Define the fixed size of the displayed portion
            self.left_display_rows = self.left_data.nrows
            self.left_display_cols = self.left_data.ncols

        # # Configure columns to expand with resizing
        self.main_ui_frame.grid_rowconfigure(
            0, weight=1
        )  # Adjust the column number if needed
        self.main_ui_frame.grid_rowconfigure(
            1, weight=1
        )  # Adjust the column number if needed
        self.main_ui_frame.grid_columnconfigure(0, weight=1)

        # Bind the click event to the canvas
        self.left_canvas.mpl_connect("button_press_event", self.on_left_canvas_click)
        self.left_canvas.mpl_connect("scroll_event", self.on_scroll)
        self.left_canvas.mpl_connect("button_release_event", self.on_left_release)
        # self.left_canvas.mpl_connect('draw_event', self.on_scroll)
        self.left_canvas.mpl_connect("motion_notify_event", self.on_canvas_motion)

    def create_right_canvas(self):
        self.right_figure, self.right_ax = plt.subplots(figsize=(6, 6))
        self.right_ax.set_title("Band Parameter Image", fontsize=8)
        self.right_ax.tick_params(axis="both", which="major", labelsize=8)
        # Set the margins and spacing to zero
        self.right_figure.tight_layout()
        self.right_figure.subplots_adjust(wspace=0, hspace=0)

        self.right_frame = tk.Frame(self.main_ui_frame)
        self.right_frame.grid(row=1, column=1, sticky="nsew")

        self.right_canvas = FigureCanvasTkAgg(
            self.right_figure, master=self.right_frame
        )
        self.right_canvas.get_tk_widget().grid(
            row=0, column=0, sticky="nsew"
        )  # Use grid instead of pack
        # self.right_canvas.get_tk_widget().grid(row=1, column=0)

        self.right_nav_toolbar_frame = tk.Frame(self.right_frame)
        self.right_nav_toolbar = NavigationToolbar2Tk(
            self.right_canvas, self.right_nav_toolbar_frame
        )
        self.right_nav_toolbar.update()
        self.right_nav_toolbar_frame.grid(row=1, column=0, sticky="nsew")

        # Configure columns to expand with resizing
        self.right_frame.grid_rowconfigure(
            0, weight=1
        )  # Adjust the column number if needed
        self.right_frame.grid_columnconfigure(0, weight=1)
        self.right_canvas.get_tk_widget().grid_columnconfigure(0, weight=1)
        self.right_canvas.get_tk_widget().grid_rowconfigure(0, weight=1)

        if self.right_data.nrows / self.right_data.ncols > 2:
            # Create variables to store the current top-left pixel coordinates
            self.top_right_row = 0
            self.top_right_col = 0

            # Define the fixed size of the displayed portion
            self.right_display_rows = int(3 * self.right_data.ncols)
            self.right_ax.set_ylim(int(3 * self.right_data.ncols), 0)
            self.right_display_cols = self.right_data.ncols
            self.right_ax.set_xlim(0, self.right_data.ncols)
        else:
            # Create variables to store the current top-left pixel coordinates
            self.top_right_row = 0
            self.top_right_col = 0

            # Define the fixed size of the displayed portion
            self.right_display_rows = self.right_data.nrows
            self.right_display_cols = self.right_data.ncols

        # Configure columns to expand with resizing
        self.main_ui_frame.grid_rowconfigure(
            0, weight=1
        )  # Adjust the column number if needed
        self.main_ui_frame.grid_rowconfigure(
            1, weight=1
        )  # Adjust the column number if needed
        self.main_ui_frame.grid_columnconfigure(1, weight=1)

        # Bind the click event to the canvas
        self.right_canvas.mpl_connect("button_press_event", self.on_left_canvas_click)
        self.right_canvas.mpl_connect("scroll_event", self.on_scroll)
        self.right_canvas.mpl_connect("button_release_event", self.on_right_release)
        # self.right_canvas.mpl_connect('button_release_event', self.on_canvas_release)
        self.right_canvas.mpl_connect("motion_notify_event", self.on_canvas_motion)

    # ----------------------------------------------------------------
    # loading the data
    # ----------------------------------------------------------------
    def load_left_data(self, file_path=None):
        if not file_path:
            file_path = filedialog.askopenfilename(
                filetypes=[
                    ("ENVI Files", "*.hdr"),
                    ("JP2 Files", "*.jp2"),
                    ("JP2 Files (Uppercase)", "*.JP2"),
                    ("TIF files", "*.tif"),
                ]
            )
            file_ext = file_path.split('.')[-1]
            print(file_ext, file_path)
        if file_path:
            # Loading a new cube invalidates anything keyed on the
            # previous source: column-lock geometry, the strip
            # footprint contour cache, and the live overlay artists.
            self._strip_footprint_cache = {}
            self._clear_strip_footprint_overlay()
            if file_ext == 'hdr' or file_ext == 'HDR':
                self.left_data = spectral.io.envi.open(file_path)
                self.ld = self.left_data.load()
                self.img_base_name = self.left_data.filename.split("/")[-1].split(".")[0]
            elif file_ext == 'tif' or file_ext == 'TIF':
                self.left_data = rio.open(file_path)
                self.left_transform = self.left_data.transform
                self.transformer = rio.transform.AffineTransformer(self.left_transform)
                self.ld = self.left_data.read()
                self.img_base_name = self.left_data.files[0].split("/")[-1].split(".")[0]
            rio_path = file_path.replace("hdr", "img")
            
            if file_ext == 'hdr' or file_ext == 'HDR':
                self.left_rio = rio.open(rio_path)
                self.left_transform = self.left_rio.transform
                self.transformer = rio.transform.AffineTransformer(self.left_transform)
            self.populate_left_wavelength_menus()
            self.create_left_canvas()
            self.display_left_data(self.default_rgb_bands)
            # elif file_ext == "jp2" or file_ext == "JP2":
            #     print('opening jp2 image')
            #     self.left_data = rio.open(file_path)
            #     self.ld = self.left_data.read().transpose(1, 2, 0)
            #     self.left_rio = rio.open(file_path)
            #     self.left_transform = self.left_rio.transform
            #     self.transformer = rio.transform.AffineTransformer(self.left_transform)
            #     self.img_base_name = file_path.split('/')[-1].split('.')[0]
            #     self.populate_left_wavelength_menus()
            #     self.create_left_canvas()
            #     self.display_left_data(self.default_rgb_bands)

            name_ = file_path.split("/")[-1].split(".")[0]
            self.root.title(f"Spectral Cube Analysis Tool: {name_}")

    def load_right_data(self, file_path=None):
        if not file_path:
            file_path = filedialog.askopenfilename(
                filetypes=[
                    ("ENVI Files", "*.hdr"),
                    ("JP2 Files", "*.jp2"),
                    ("JP2 Files (Uppercase)", "*.JP2"),
                ]
            )
        if file_path:
            self.right_data = spectral.io.envi.open(file_path)
            self.rd = self.right_data.load()
            self.populate_right_wavelength_menus()
            self.create_right_canvas()
            self.display_right_data(self.default_parameter_bands)

    def populate_left_wavelength_menus(self):
        self.bands_inverted = False
        if self.left_data is not None:
            try:
                if self.left_data.bands.centers is not None:
                    # Check if parameter image or spectral cube
                    if str(self.left_data.bands.centers[1]).isalpha():
                        # This is a band parameter image
                        self.left_is_parameter = True
                        wavelengths = self.left_data.bands.centers
                        self.left_wvl = self.left_data.bands.centers
                    else:
                        # This is a spectral cube
                        if (
                            self.left_data.bands.centers[1]
                            - self.left_data.bands.centers[-1]
                            > 0
                        ):
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.left_data.bands.centers
                            self.left_wvl = sorted(
                                [float(i) for i in self.left_data.bands.centers]
                            )
                        else:
                            wavelengths = self.left_data.bands.centers
                            self.left_wvl = [
                                float(i) for i in self.left_data.bands.centers
                            ]

                elif self.left_data.metadata["band names"] is not None:
                    if str(self.left_data.metadata["band names"][1]).isalpha():
                        # This is a band parameter image
                        self.left_is_parameter = True
                        wavelengths = self.left_data.metadata["band names"]
                        self.left_wvl = self.left_data.metadata["band names"]
                    else:
                        # This is a spectral cube
                        wavelengths = (
                            np.array(self.left_data.metadata["band names"])
                            .astype(float)
                            .tolist()
                        )
                        self.left_data.metadata["band names"] = wavelengths
                        if (
                            self.left_data.metadata["band names"][1]
                            - self.left_data.metadata["band names"][-1]
                            > 0
                        ):
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.left_data.metadata["band names"]
                            self.left_wvl = sorted(
                                [
                                    float(i)
                                    for i in self.left_data.metadata["band names"]
                                ]
                            )
                        else:
                            # This is a band parameter image
                            wavelengths = self.left_data.metadata["band names"]
                            self.left_wvl = self.left_data.metadata["band names"]
            except:
                wavelengths = list(range(1, self.left_data.shape[2] + 1))
                self.left_wvl = wavelengths
                messagebox.showerror(
                    "Error",
                    "Unable to load wavelength information. Setting wavelengths to default integer values",
                )
            

            # set the wavelength values
            self.left_red_band_menu["values"] = wavelengths
            self.left_green_band_menu["values"] = wavelengths
            self.left_blue_band_menu["values"] = wavelengths

            # Find default band indices based on target wavelengths or band names
            if self.left_is_parameter:
                default_indices = self.get_default_parameter_indices(wavelengths)
            else:
                default_indices = self.get_default_spectral_indices(wavelengths)

            self.left_red_band_menu.set(wavelengths[default_indices[0]])
            self.left_green_band_menu.set(wavelengths[default_indices[1]])
            self.left_blue_band_menu.set(wavelengths[default_indices[2]])

            # Calculate and set stretch limits for each channel
            red_index, green_index, blue_index = default_indices
            method = self.left_stretch_type_var.get()

            red_min_stretch, red_max_stretch = self.calculate_stretch_limits(
                self.left_data[:, :, red_index], method
            )
            green_min_stretch, green_max_stretch = self.calculate_stretch_limits(
                self.left_data[:, :, green_index], method
            )
            blue_min_stretch, blue_max_stretch = self.calculate_stretch_limits(
                self.left_data[:, :, blue_index], method
            )

            self.left_red_min_stretch_var.set(red_min_stretch)
            self.left_red_max_stretch_var.set(red_max_stretch)
            self.left_green_min_stretch_var.set(green_min_stretch)
            self.left_green_max_stretch_var.set(green_max_stretch)
            self.left_blue_min_stretch_var.set(blue_min_stretch)
            self.left_blue_max_stretch_var.set(blue_max_stretch)

    def populate_right_wavelength_menus(self):
        self.bands_inverted = False
        if self.right_data is not None:
            try:
                if self.right_data.bands.centers is not None:
                    # check if parameter image or spectral cube
                    if str(self.right_data.bands.centers[1]).isalpha():
                        # this is a band parameter image
                        self.right_is_parameter = True
                        wavelengths = self.right_data.bands.centers
                        self.right_wvl = self.right_data.bands.centers
                    else:
                        # this is a spectral cube
                        if (
                            self.right_data.bands.centers[1]
                            - self.right_data.bands.centers[-1]
                            > 0
                        ):
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.right_data.bands.centers
                            self.right_wvl = sorted(
                                [float(i) for i in self.right_data.bands.centers]
                            )
                        else:
                            wavelengths = self.right_data.bands.centers
                            self.right_wvl = [
                                float(i) for i in self.right_data.bands.centers
                            ]

                elif (
                    self.right_data.bands.centers is None
                    and self.right_data.metadata["band names"] is not None
                ):
                    if any(
                        char.isalpha()
                        for char in str(self.right_data.metadata["band names"][1])
                    ):
                        # this is a band parameter image
                        self.right_is_parameter = True
                        wavelengths = self.right_data.metadata["band names"]
                        self.right_wvl = self.right_data.metadata["band names"]
                    else:
                        # this is a spectral cube
                        if (
                            self.right_data.metadata["band names"][1]
                            - self.right_data.metadata["band names"][-1]
                            > 0
                        ):
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.right_data.metadata["band names"]
                            self.right_wvl = sorted(
                                [
                                    float(i)
                                    for i in self.right_data.metadata["band names"]
                                ]
                            )
                        else:
                            # this is a band parameter image
                            wavelengths = self.right_data.metadata["band names"]
                            self.right_wvl = self.right_data.metadata["band names"]
            except:
                wavelengths = list(range(1, self.right_data.shape[2] + 1))
                self.right_wvl = wavelengths
                messagebox.showerror("Error", "Unable to load wavelength information.")

            self.right_red_band_menu["values"] = wavelengths
            self.right_green_band_menu["values"] = wavelengths
            self.right_blue_band_menu["values"] = wavelengths

            # Find default band indices based on target band names or wavelengths
            if self.right_is_parameter:
                default_indices = self.get_default_parameter_indices(wavelengths)
            else:
                default_indices = self.get_default_spectral_indices(wavelengths)

            self.right_red_band_menu.set(wavelengths[default_indices[0]])
            self.right_green_band_menu.set(wavelengths[default_indices[1]])
            self.right_blue_band_menu.set(wavelengths[default_indices[2]])

            # Calculate and set stretch limits for each channel in the "right_data"
            red_index, green_index, blue_index = default_indices
            method = self.right_stretch_type_var.get()

            right_red_min_stretch, right_red_max_stretch = self.calculate_stretch_limits(
                self.right_data[:, :, red_index], method
            )
            right_green_min_stretch, right_green_max_stretch = self.calculate_stretch_limits(
                self.right_data[:, :, green_index], method
            )
            right_blue_min_stretch, right_blue_max_stretch = self.calculate_stretch_limits(
                self.right_data[:, :, blue_index], method
            )

            self.right_red_min_stretch_var.set(right_red_min_stretch)
            self.right_red_max_stretch_var.set(right_red_max_stretch)
            self.right_green_min_stretch_var.set(right_green_min_stretch)
            self.right_green_max_stretch_var.set(right_green_max_stretch)
            self.right_blue_min_stretch_var.set(right_blue_min_stretch)
            self.right_blue_max_stretch_var.set(right_blue_max_stretch)

    # ----------------------------------------------------------------
    # displaying the data
    # ----------------------------------------------------------------
    def update_left_display(self):
        try:
            red_band = float(self.left_red_band_var.get())
            green_band = float(self.left_green_band_var.get())
            blue_band = float(self.left_blue_band_var.get())
            red_index = self.left_wvl.index(red_band)
            green_index = self.left_wvl.index(green_band)
            blue_index = self.left_wvl.index(blue_band)
            self.left_band_indices = [red_index, green_index, blue_index]
            self.display_left_data(self.left_band_indices)

            if (
                self.left_hist_window is not None
                and self.left_hist_window.winfo_exists()
            ):
                self.plot_left_histograms()

        except ValueError:
            messagebox.showerror(
                f"Error",
                "Invalid wavelength. Please enter valid wavelengths.\n{ValueError}",
            )

    def update_right_display(self):
        try:
            # check if it's a spectral cube or band parameter image
            if self.right_is_parameter:
                red_band = self.right_red_band_var.get()
                green_band = self.right_green_band_var.get()
                blue_band = self.right_blue_band_var.get()
            else:
                red_band = float(self.right_red_band_var.get())
                green_band = float(self.right_green_band_var.get())
                blue_band = float(self.right_blue_band_var.get())
            red_index = self.right_wvl.index(red_band)
            green_index = self.right_wvl.index(green_band)
            blue_index = self.right_wvl.index(blue_band)
            self.right_band_indices = [red_index, green_index, blue_index]

            self.display_right_data(self.right_band_indices)

            if (
                self.right_hist_window is not None
                and self.right_hist_window.winfo_exists()
            ):
                self.plot_right_histograms()

        except ValueError:
            messagebox.showerror(
                f"Error",
                "Invalid wavelength. Please enter valid wavelengths.\n{ValueError}",
            )

    def display_left_data(self, band_indices):

        left_red_stretch = (
            (self.left_red_min_stretch_var.get(), self.left_red_max_stretch_var.get())
            if hasattr(self, "left_red_min_stretch_var")
            else None
        )
        left_green_stretch = (
            (
                self.left_green_min_stretch_var.get(),
                self.left_green_max_stretch_var.get(),
            )
            if hasattr(self, "left_green_min_stretch_var")
            else None
        )
        left_blue_stretch = (
            (self.left_blue_min_stretch_var.get(), self.left_blue_max_stretch_var.get())
            if hasattr(self, "left_blue_min_stretch_var")
            else None
        )

        # print(self.top_left_row, self.left_display_rows, self.top_left_col, self.left_display_cols)
        # left_rgb_image = self.left_data[self.top_left_row:self.top_left_row + self.left_display_rows,
        #                                 self.top_left_col:self.top_left_col + self.left_display_cols,
        #                                 [band_indices[0], band_indices[1], band_indices[2]]]

        left_rgb_image = self.left_data[
            :, :, [band_indices[0], band_indices[1], band_indices[2]]
        ]

        method = self.left_stretch_type_var.get() if hasattr(self, 'left_stretch_type_var') else "linear"
        if left_red_stretch is not None:
            left_rgb_image[:, :, 0] = self.stretch_band(
                left_rgb_image[:, :, 0], left_red_stretch, method
            )
        if left_green_stretch is not None:
            left_rgb_image[:, :, 1] = self.stretch_band(
                left_rgb_image[:, :, 1], left_green_stretch, method
            )
        if left_blue_stretch is not None:
            left_rgb_image[:, :, 2] = self.stretch_band(
                left_rgb_image[:, :, 2], left_blue_stretch, method
            )

        self.left_ax.imshow(np.array(left_rgb_image))
        self.left_canvas.draw()
        # show(np.transpose(left_rgb_image, (2,0,1)), ax=self.left_ax, transform=self.left_transform)

    def display_right_data(self, band_indices):

        right_red_stretch = (
            (self.right_red_min_stretch_var.get(), self.right_red_max_stretch_var.get())
            if hasattr(self, "right_red_min_stretch_var")
            else None
        )
        right_green_stretch = (
            (
                self.right_green_min_stretch_var.get(),
                self.right_green_max_stretch_var.get(),
            )
            if hasattr(self, "right_green_min_stretch_var")
            else None
        )
        right_blue_stretch = (
            (
                self.right_blue_min_stretch_var.get(),
                self.right_blue_max_stretch_var.get(),
            )
            if hasattr(self, "right_blue_min_stretch_var")
            else None
        )

        # right_rgb_image = self.right_data[self.top_right_row:self.top_right_row + self.right_display_rows,
        #                                 self.top_right_col:self.top_right_col + self.right_display_cols,
        #                                 [band_indices[0], band_indices[1], band_indices[2]]]

        right_rgb_image = self.right_data[
            :, :, [band_indices[0], band_indices[1], band_indices[2]]
        ]

        method = self.right_stretch_type_var.get() if hasattr(self, 'right_stretch_type_var') else "linear"
        if right_red_stretch is not None:
            right_rgb_image[:, :, 0] = self.stretch_band(
                right_rgb_image[:, :, 0], right_red_stretch, method
            )
        if right_green_stretch is not None:
            right_rgb_image[:, :, 1] = self.stretch_band(
                right_rgb_image[:, :, 1], right_green_stretch, method
            )
        if right_blue_stretch is not None:
            right_rgb_image[:, :, 2] = self.stretch_band(
                right_rgb_image[:, :, 2], right_blue_stretch, method
            )

        self.right_ax.imshow(np.array(right_rgb_image))
        self.right_canvas.draw()

    def reset_display_extent(self):
        self.left_ax.set_ylim(self.left_display_rows, 0)
        self.right_ax.set_ylim(self.right_display_rows, 0)
        self.left_ax.set_xlim(0, self.left_display_cols)
        self.right_ax.set_xlim(0, self.right_display_cols)
        self.left_canvas.draw()
        self.right_canvas.draw()

    def apply_new_left_bands(self):
        red_band = float(self.left_red_band_var.get())
        green_band = float(self.left_green_band_var.get())
        blue_band = float(self.left_blue_band_var.get())
        red_index = self.left_wvl.index(red_band)
        green_index = self.left_wvl.index(green_band)
        blue_index = self.left_wvl.index(blue_band)

        method = self.left_stretch_type_var.get()
        r_min, r_max = self.calculate_stretch_limits(self.left_data[:, :, red_index], method)
        g_min, g_max = self.calculate_stretch_limits(self.left_data[:, :, green_index], method)
        b_min, b_max = self.calculate_stretch_limits(self.left_data[:, :, blue_index], method)

        self.left_red_min_stretch_var.set(r_min)
        self.left_red_max_stretch_var.set(r_max)
        self.left_green_min_stretch_var.set(g_min)
        self.left_green_max_stretch_var.set(g_max)
        self.left_blue_min_stretch_var.set(b_min)
        self.left_blue_max_stretch_var.set(b_max)

        self.update_left_display()

    def apply_new_right_bands(self):
        # check if it's a spectral cube or band parameter image
        if self.right_is_parameter:
            red_band = self.right_red_band_var.get()
            green_band = self.right_green_band_var.get()
            blue_band = self.right_blue_band_var.get()
        else:
            red_band = float(self.right_red_band_var.get())
            green_band = float(self.right_green_band_var.get())
            blue_band = float(self.right_blue_band_var.get())
        red_index = self.right_wvl.index(red_band)
        green_index = self.right_wvl.index(green_band)
        blue_index = self.right_wvl.index(blue_band)

        method = self.right_stretch_type_var.get()
        r_min, r_max = self.calculate_stretch_limits(self.right_data[:, :, red_index], method)
        g_min, g_max = self.calculate_stretch_limits(self.right_data[:, :, green_index], method)
        b_min, b_max = self.calculate_stretch_limits(self.right_data[:, :, blue_index], method)

        self.right_red_min_stretch_var.set(r_min)
        self.right_red_max_stretch_var.set(r_max)
        self.right_green_min_stretch_var.set(g_min)
        self.right_green_max_stretch_var.set(g_max)
        self.right_blue_min_stretch_var.set(b_min)
        self.right_blue_max_stretch_var.set(b_max)

        self.update_right_display()

    # ----------------------------------------------------------------
    # histograms and stretching
    # ----------------------------------------------------------------
    def stretch_band(self, band, stretch_range, method="linear"):
        if method == "Hist. Equalization":
            return self.histogram_equalize_band(band)
        min_val, max_val = stretch_range
        if min_val >= max_val:
            # Prevent division by zero or inverted stretch
            max_val = min_val + 1e-10
        stretched_band = (band - min_val) / (max_val - min_val)
        stretched_band = np.clip(stretched_band, 0, 1)
        return stretched_band

    def histogram_equalize_band(self, band):
        # Flatten and get valid pixels (exclude zeros, NaNs, and CRISM fill value 65535)
        flat = band.flatten()
        valid_mask = (np.isfinite(flat) & (flat != 0) &
                      (flat != 65535.0) & (flat != 65535))
        valid_data = flat[valid_mask]
        if len(valid_data) == 0:
            return np.zeros_like(band)
        # Build CDF from valid data
        sorted_data = np.sort(valid_data)
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
        # Map all pixels using the CDF
        equalized = np.interp(band, sorted_data, cdf)
        # Set invalid pixels to 0
        invalid_mask = (~np.isfinite(band) | (band == 0) |
                        (band == 65535.0) | (band == 65535))
        equalized[invalid_mask] = 0
        return equalized

    def calculate_stretch_limits(self, channel_data, method="Percentile"):
        # Filter out zeros (padding), NaNs, and CRISM fill value (65535)
        valid_mask = (np.isfinite(channel_data) &
                      (channel_data != 0) &
                      (channel_data != 65535.0) &
                      (channel_data != 65535))
        valid_data = channel_data[valid_mask]

        if len(valid_data) == 0:
            return 0.0, 1.0

        if method == "Percentile":
            min_stretch = np.percentile(valid_data, 2)
            max_stretch = np.percentile(valid_data, 98)
        elif method == "Sigma Clip":
            # Iterative sigma clipping: remove 3-sigma outliers, then use mean +/- 2*std
            clipped = valid_data.copy()
            for _ in range(3):
                mean_val = np.mean(clipped)
                std_val = np.std(clipped)
                if std_val == 0:
                    break
                mask = np.abs(clipped - mean_val) < 3 * std_val
                if not np.any(mask):
                    break
                clipped = clipped[mask]
            mean_val = np.mean(clipped)
            std_val = np.std(clipped)
            min_stretch = mean_val - 2 * std_val
            max_stretch = mean_val + 2 * std_val
        else:
            # Fallback to percentile
            min_stretch = np.percentile(valid_data, 2)
            max_stretch = np.percentile(valid_data, 98)

        # Ensure min < max
        if min_stretch >= max_stretch:
            max_stretch = min_stretch + 1e-10

        return float(min_stretch), float(max_stretch)

    def auto_stretch_left(self):
        if not hasattr(self, 'left_data') or self.left_data is None:
            return
        method = self.left_stretch_type_var.get()
        if method == "Hist. Equalization":
            # For equalization, just update display directly
            self.update_left_display()
            return
        red_index = self.left_wvl.index(float(self.left_red_band_var.get()))
        green_index = self.left_wvl.index(float(self.left_green_band_var.get()))
        blue_index = self.left_wvl.index(float(self.left_blue_band_var.get()))

        r_min, r_max = self.calculate_stretch_limits(self.left_data[:, :, red_index], method)
        g_min, g_max = self.calculate_stretch_limits(self.left_data[:, :, green_index], method)
        b_min, b_max = self.calculate_stretch_limits(self.left_data[:, :, blue_index], method)

        self.left_red_min_stretch_var.set(r_min)
        self.left_red_max_stretch_var.set(r_max)
        self.left_green_min_stretch_var.set(g_min)
        self.left_green_max_stretch_var.set(g_max)
        self.left_blue_min_stretch_var.set(b_min)
        self.left_blue_max_stretch_var.set(b_max)
        self.update_left_display()

    def auto_stretch_right(self):
        if not hasattr(self, 'right_data') or self.right_data is None:
            return
        method = self.right_stretch_type_var.get()
        if method == "Hist. Equalization":
            self.update_right_display()
            return
        if self.right_is_parameter:
            red_band = self.right_red_band_var.get()
            green_band = self.right_green_band_var.get()
            blue_band = self.right_blue_band_var.get()
            red_index = list(self.right_data.bands.centers).index(red_band)
            green_index = list(self.right_data.bands.centers).index(green_band)
            blue_index = list(self.right_data.bands.centers).index(blue_band)
        else:
            red_index = self.right_wvl.index(float(self.right_red_band_var.get()))
            green_index = self.right_wvl.index(float(self.right_green_band_var.get()))
            blue_index = self.right_wvl.index(float(self.right_blue_band_var.get()))

        r_min, r_max = self.calculate_stretch_limits(self.right_data[:, :, red_index], method)
        g_min, g_max = self.calculate_stretch_limits(self.right_data[:, :, green_index], method)
        b_min, b_max = self.calculate_stretch_limits(self.right_data[:, :, blue_index], method)

        self.right_red_min_stretch_var.set(r_min)
        self.right_red_max_stretch_var.set(r_max)
        self.right_green_min_stretch_var.set(g_min)
        self.right_green_max_stretch_var.set(g_max)
        self.right_blue_min_stretch_var.set(b_min)
        self.right_blue_max_stretch_var.set(b_max)
        self.update_right_display()

    def find_closest_wavelength_index(self, wavelengths, target_wvl):
        """Find the index of the wavelength closest to target_wvl."""
        wvl_array = np.array([float(w) for w in wavelengths])
        idx = np.argmin(np.abs(wvl_array - target_wvl))
        return int(idx)

    def find_band_name_index(self, band_names, target_name):
        """Find the index of a band by name (case-insensitive partial match)."""
        target_upper = target_name.upper()
        for i, name in enumerate(band_names):
            if target_upper in str(name).upper():
                return i
        return None

    def get_default_spectral_indices(self, wavelengths):
        """Get indices for default CRISM RGB wavelengths (2330, 1500, 774 nm)."""
        r_idx = self.find_closest_wavelength_index(wavelengths, self.default_rgb_wavelengths[0])
        g_idx = self.find_closest_wavelength_index(wavelengths, self.default_rgb_wavelengths[1])
        b_idx = self.find_closest_wavelength_index(wavelengths, self.default_rgb_wavelengths[2])
        return [r_idx, g_idx, b_idx]

    def get_default_parameter_indices(self, band_names):
        """Get indices for default CRISM parameter bands (BD1300, RPEAK1, LCPINDEX2)."""
        indices = []
        for target in self.default_parameter_names:
            idx = self.find_band_name_index(band_names, target)
            if idx is not None:
                indices.append(idx)
            else:
                # Fallback: use position in default list
                fallback_idx = len(indices)
                if fallback_idx < len(band_names):
                    indices.append(fallback_idx)
                else:
                    indices.append(0)
        return indices

    def plot_left_histograms(self):
        if hasattr(self, "left_data"):
            left_red_band = float(self.left_red_band_var.get())
            left_green_band = float(self.left_green_band_var.get())
            left_blue_band = float(self.left_blue_band_var.get())
            left_red_index = self.left_wvl.index(left_red_band)
            left_green_index = self.left_wvl.index(left_green_band)
            left_blue_index = self.left_wvl.index(left_blue_band)

            left_red_channel = np.where(
                self.left_data[:, :, left_red_index] > 1,
                np.nan,
                self.left_data[:, :, left_red_index],
            )
            left_green_channel = np.where(
                self.left_data[:, :, left_green_index] > 1,
                np.nan,
                self.left_data[:, :, left_green_index],
            )
            left_blue_channel = np.where(
                self.left_data[:, :, left_blue_index] > 1,
                np.nan,
                self.left_data[:, :, left_blue_index],
            )

            # Calculate min and max values across all three channels
            left_min_value = min(
                np.nanmin(left_red_channel),
                np.nanmin(left_green_channel),
                np.nanmin(left_blue_channel),
            )
            left_max_value = max(
                np.nanmax(left_red_channel),
                np.nanmax(left_green_channel),
                np.nanmax(left_blue_channel),
            )

            if (
                self.left_hist_window is None
                or not self.left_hist_window.winfo_exists()
            ):
                hist_args = (
                    left_min_value,
                    left_max_value,
                    left_red_band,
                    left_green_band,
                    left_blue_band,
                )
                self.create_left_histogram(hist_args)

            if self.reset_left_hist:
                red_range = (left_min_value, left_max_value)
                green_range = (left_min_value, left_max_value)
                blue_range = (left_min_value, left_max_value)
                self.reset_left_hist = False
            else:
                red_range = (
                    self.left_red_min_stretch_var.get()
                    - abs(self.left_red_min_stretch_var.get() * 0.2),
                    self.left_red_max_stretch_var.get()
                    + abs(self.left_red_max_stretch_var.get() * 0.2),
                )
                green_range = (
                    self.left_green_min_stretch_var.get()
                    - abs(self.left_green_min_stretch_var.get() * 0.2),
                    self.left_green_max_stretch_var.get()
                    + abs(self.left_green_max_stretch_var.get() * 0.2),
                )
                blue_range = (
                    self.left_blue_min_stretch_var.get()
                    - abs(self.left_blue_min_stretch_var.get() * 0.2),
                    self.left_blue_max_stretch_var.get()
                    + abs(self.left_blue_max_stretch_var.get() * 0.2),
                )

            # Update the data in the existing subplots
            self.left_hist_axes[0].cla()
            self.left_hist_axes[0].hist(
                left_red_channel.ravel(),
                bins=256,
                color="red",
                alpha=0.7,
                range=red_range,
            )
            self.left_hist_axes[0].axvline(
                self.left_red_min_stretch_var.get(),
                color="red",
                linestyle="--",
                label="Min Stretch",
            )
            self.left_hist_axes[0].axvline(
                self.left_red_max_stretch_var.get(),
                color="red",
                linestyle="--",
                label="Max Stretch",
            )
            self.left_hist_axes[0].set_title(
                f"Left Red Channel Histogram: {left_red_band}"
            )
            # self.left_hist_axes[0].legend()

            self.left_hist_axes[1].cla()
            self.left_hist_axes[1].hist(
                left_green_channel.ravel(),
                bins=256,
                color="green",
                alpha=0.7,
                range=green_range,
            )
            self.left_hist_axes[1].axvline(
                self.left_green_min_stretch_var.get(),
                color="green",
                linestyle="--",
                label="Min Stretch",
            )
            self.left_hist_axes[1].axvline(
                self.left_green_max_stretch_var.get(),
                color="green",
                linestyle="--",
                label="Max Stretch",
            )
            self.left_hist_axes[1].set_title(
                f"Left Green Channel Histogram: {left_green_band}"
            )
            # self.left_hist_axes[1].legend()

            self.left_hist_axes[2].cla()
            self.left_hist_axes[2].hist(
                left_blue_channel.ravel(),
                bins=256,
                color="blue",
                alpha=0.7,
                range=blue_range,
            )
            self.left_hist_axes[2].axvline(
                self.left_blue_min_stretch_var.get(),
                color="blue",
                linestyle="--",
                label="Min Stretch",
            )
            self.left_hist_axes[2].axvline(
                self.left_blue_max_stretch_var.get(),
                color="blue",
                linestyle="--",
                label="Max Stretch",
            )
            self.left_hist_axes[2].set_title(
                f"Left Blue Channel Histogram: {left_blue_band}"
            )
            # self.left_hist_axes[2].legend()

            # Update the existing figure
            self.left_hist_figure.canvas.draw()
            self.left_hist_figure.canvas.flush_events()

            # Create draggable span selectors
            self.create_left_histogram_span_selectors()

        else:
            messagebox.showwarning(
                "Warning",
                "No left frame data loaded. Load data into the left frame first.",
            )

    def plot_right_histograms(self):
        if hasattr(self, "right_data"):
            # check if it's a spectral cube or band parameter image
            if self.right_is_parameter:
                right_red_band = self.right_red_band_var.get()
                right_green_band = self.right_green_band_var.get()
                right_blue_band = self.right_blue_band_var.get()
            else:
                right_red_band = float(self.right_red_band_var.get())
                right_green_band = float(self.right_green_band_var.get())
                right_blue_band = float(self.right_blue_band_var.get())

            right_red_index = self.right_wvl.index(right_red_band)
            right_green_index = self.right_wvl.index(right_green_band)
            right_blue_index = self.right_wvl.index(right_blue_band)

            # could update this to plot the histogram of the zoom section
            # right_red_channel = np.where(self.right_data[:, :, right_red_index] > 1, np.nan, self.right_data[:, :, right_red_index])
            # right_green_channel = np.where(self.right_data[:, :, right_green_index] > 1, np.nan, self.right_data[:, :, right_green_index])
            # right_blue_channel = np.where(self.right_data[:, :, right_blue_index] > 1, np.nan, self.right_data[:, :, right_blue_index])
            right_red_channel = self.right_data[:, :, right_red_index]
            right_green_channel = self.right_data[:, :, right_green_index]
            right_blue_channel = self.right_data[:, :, right_blue_index]

            # Calculate min and max values across all three channels
            # right_min_value = min(np.nanmin(right_red_channel), np.nanmin(right_green_channel), np.nanmin(right_blue_channel))
            # right_max_value = max(np.nanmax(right_red_channel), np.nanmax(right_green_channel), np.nanmax(right_blue_channel))

            # calculate the mean and standard deviation of all three channels
            right_red_mean = np.nanmedian(right_red_channel)
            right_green_mean = np.nanmedian(right_green_channel)
            right_blue_mean = np.nanmedian(right_blue_channel)
            right_red_std = np.nanquantile(right_red_channel, 0.6)
            right_green_std = np.nanquantile(right_green_channel, 0.6)
            right_blue_std = np.nanquantile(right_blue_channel, 0.6)

            # set the right_min_value and right_max_value to be -3*std and 3*std from the mean
            right_min_value = min(
                right_red_mean - right_red_std,
                right_green_mean - right_green_std,
                right_blue_mean - right_blue_std,
            )
            right_max_value = max(
                right_red_mean + right_red_std,
                right_green_mean + right_green_std,
                right_blue_mean + right_blue_std,
            )

            if (
                self.right_hist_window is None
                or not self.right_hist_window.winfo_exists()
            ):
                hist_args = (
                    right_min_value,
                    right_max_value,
                    right_red_band,
                    right_green_band,
                    right_blue_band,
                )
                self.create_right_histogram(hist_args)

            if self.reset_right_hist:
                red_range = (right_min_value, right_max_value)
                green_range = (right_min_value, right_max_value)
                blue_range = (right_min_value, right_max_value)
                self.reset_right_hist = False
            else:
                red_range = (
                    self.right_red_min_stretch_var.get()
                    - abs(self.right_red_min_stretch_var.get() * 0.2),
                    self.right_red_max_stretch_var.get()
                    + abs(self.right_red_max_stretch_var.get() * 0.2),
                )
                green_range = (
                    self.right_green_min_stretch_var.get()
                    - abs(self.right_green_min_stretch_var.get() * 0.2),
                    self.right_green_max_stretch_var.get()
                    + abs(self.right_green_max_stretch_var.get() * 0.2),
                )
                blue_range = (
                    self.right_blue_min_stretch_var.get()
                    - abs(self.right_blue_min_stretch_var.get() * 0.2),
                    self.right_blue_max_stretch_var.get()
                    + abs(self.right_blue_max_stretch_var.get() * 0.2),
                )

            # Update the data in the existing subplots
            self.right_hist_axes[0].cla()
            self.right_hist_axes[0].hist(
                right_red_channel.ravel(),
                bins=256,
                color="red",
                alpha=0.7,
                range=red_range,
            )
            self.right_hist_axes[0].axvline(
                self.right_red_min_stretch_var.get(),
                color="red",
                linestyle="--",
                label="Min Stretch",
            )
            self.right_hist_axes[0].axvline(
                self.right_red_max_stretch_var.get(),
                color="red",
                linestyle="--",
                label="Max Stretch",
            )
            self.right_hist_axes[0].set_title(
                f"Right Red Channel Histogram: {right_red_band}"
            )
            # self.right_hist_axes[0].legend()

            self.right_hist_axes[1].cla()
            self.right_hist_axes[1].hist(
                right_green_channel.ravel(),
                bins=256,
                color="green",
                alpha=0.7,
                range=green_range,
            )
            self.right_hist_axes[1].axvline(
                self.right_green_min_stretch_var.get(),
                color="green",
                linestyle="--",
                label="Min Stretch",
            )
            self.right_hist_axes[1].axvline(
                self.right_green_max_stretch_var.get(),
                color="green",
                linestyle="--",
                label="Max Stretch",
            )
            self.right_hist_axes[1].set_title(
                f"Right Green Channel Histogram: {right_green_band}"
            )
            # self.right_hist_axes[1].legend()

            self.right_hist_axes[2].cla()
            self.right_hist_axes[2].hist(
                right_blue_channel.ravel(),
                bins=256,
                color="blue",
                alpha=0.7,
                range=blue_range,
            )
            self.right_hist_axes[2].axvline(
                self.right_blue_min_stretch_var.get(),
                color="blue",
                linestyle="--",
                label="Min Stretch",
            )
            self.right_hist_axes[2].axvline(
                self.right_blue_max_stretch_var.get(),
                color="blue",
                linestyle="--",
                label="Max Stretch",
            )
            self.right_hist_axes[2].set_title(
                f"Right Blue Channel Histogram: {right_blue_band}"
            )
            # self.right_hist_axes[2].legend()

            # Update the existing figure
            self.right_hist_figure.canvas.draw()
            self.right_hist_figure.canvas.flush_events()

            # Create draggable span selectors
            self.create_right_histogram_span_selectors()

        else:
            messagebox.showwarning(
                "Warning",
                "No right frame data loaded. Load data into the right frame first.",
            )

    def reset_left_hist_button(self):
        self.reset_left_hist = True
        self.plot_left_histograms()

    def create_left_histogram(self, hist_args):
        (
            left_min_value,
            left_max_value,
            left_red_band,
            left_green_band,
            left_blue_band,
        ) = hist_args
        self.left_hist_window = tk.Toplevel(self.root)
        self.left_hist_window.title("Hyperspectral Image Histogram")
        # self.left_hist_window.geometry("600x800")

        # Create a frame to hold UI elements with a fixed size
        left_hist_ui_frame = tk.Frame(self.left_hist_window)
        left_hist_ui_frame.pack(fill=tk.X)

        # Create a new histogram plot window if it doesn't exist
        self.left_hist_figure, self.left_hist_axes = plt.subplots(3, 1, figsize=(5, 6))

        self.left_hist_canvas = FigureCanvasTkAgg(
            self.left_hist_figure, master=self.left_hist_window
        )
        self.left_hist_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset the histogram
        self.reset_left_hist_x_axis_button = tk.Button(
            left_hist_ui_frame,
            text="Reset Histogram",
            command=self.reset_left_hist_button,
        )
        self.reset_left_hist_x_axis_button.pack(side=tk.RIGHT)

        # Create the subplots within the figure
        self.left_hist_axes[0].hist(
            [], bins=256, color="red", alpha=0.7, range=(left_min_value, left_max_value)
        )
        self.left_hist_axes[0].set_title(f"Left Red Channel Histogram: {left_red_band}")
        self.left_hist_axes[0].legend()

        self.left_hist_axes[1].hist(
            [],
            bins=256,
            color="green",
            alpha=0.7,
            range=(left_min_value, left_max_value),
        )
        self.left_hist_axes[1].set_title(
            f"Left Green Channel Histogram: {left_green_band}"
        )
        self.left_hist_axes[1].legend()

        self.left_hist_axes[2].hist(
            [],
            bins=256,
            color="blue",
            alpha=0.7,
            range=(left_min_value, left_max_value),
        )
        self.left_hist_axes[2].set_title(
            f"Left Blue Channel Histogram: {left_blue_band}"
        )
        self.left_hist_axes[2].legend()

        # Add spacing between subplots
        self.left_hist_figure.subplots_adjust(hspace=0.8)

        # Create a toolbar for the spectral plot
        left_hist_toolbar = NavigationToolbar2Tk(
            self.left_hist_canvas, self.left_hist_window
        )
        left_hist_toolbar.update()
        self.left_hist_canvas.get_tk_widget().pack(
            side=tk.TOP, fill=tk.BOTH, expand=True
        )

    def reset_right_hist_button(self):
        self.reset_right_hist = True
        self.plot_right_histograms()

    def create_right_histogram(self, hist_args):
        (
            right_min_value,
            right_max_value,
            right_red_band,
            right_green_band,
            right_blue_band,
        ) = hist_args
        self.right_hist_window = tk.Toplevel(self.root)
        self.right_hist_window.title("Band Parameter Histogram")
        # self.right_hist_window.geometry("600x800")

        # Create a frame to hold UI elements with a fixed size
        right_hist_ui_frame = tk.Frame(self.right_hist_window)
        right_hist_ui_frame.pack(fill=tk.X)

        # Create a new histogram plot window if it doesn't exist
        self.right_hist_figure, self.right_hist_axes = plt.subplots(
            3, 1, figsize=(5, 6)
        )

        self.right_hist_canvas = FigureCanvasTkAgg(
            self.right_hist_figure, master=self.right_hist_window
        )
        self.right_hist_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset the histogram
        self.reset_right_hist_x_axis_button = tk.Button(
            right_hist_ui_frame,
            text="Reset Histogram",
            command=self.reset_right_hist_button,
        )
        self.reset_right_hist_x_axis_button.pack(side=tk.RIGHT)

        # Create the subplots within the figure
        self.right_hist_axes[0].hist(
            [],
            bins=256,
            color="red",
            alpha=0.7,
            range=(right_min_value, right_max_value),
        )
        self.right_hist_axes[0].set_title(
            f"Right Red Channel Histogram: {right_red_band}", fontsize=4
        )
        self.right_hist_axes[0].legend()

        self.right_hist_axes[1].hist(
            [],
            bins=256,
            color="green",
            alpha=0.7,
            range=(right_min_value, right_max_value),
        )
        self.right_hist_axes[1].set_title(
            f"Right Green Channel Histogram: {right_green_band}", fontsize=4
        )
        self.right_hist_axes[1].legend()

        self.right_hist_axes[2].hist(
            [],
            bins=256,
            color="blue",
            alpha=0.7,
            range=(right_min_value, right_max_value),
        )
        self.right_hist_axes[2].set_title(
            f"Right Blue Channel Histogram: {right_blue_band}", fontsize=4
        )
        self.right_hist_axes[2].legend()

        # Add spacing between subplots
        self.right_hist_figure.subplots_adjust(hspace=0.8)

        # Create a toolbar for the spectral plot
        right_hist_toolbar = NavigationToolbar2Tk(
            self.right_hist_canvas, self.right_hist_window
        )
        right_hist_toolbar.update()
        self.right_hist_canvas.get_tk_widget().pack(
            side=tk.TOP, fill=tk.BOTH, expand=True
        )

    def update_left_red_stretch(self, val_min, val_max):
        self.left_red_min_stretch_var.set(val_min)
        self.left_red_max_stretch_var.set(val_max)
        self.plot_left_histograms()
        self.update_left_display()

    def update_left_green_stretch(self, val_min, val_max):
        self.left_green_min_stretch_var.set(val_min)
        self.left_green_max_stretch_var.set(val_max)
        self.plot_left_histograms()
        self.update_left_display()

    def update_left_blue_stretch(self, val_min, val_max):
        self.left_blue_min_stretch_var.set(val_min)
        self.left_blue_max_stretch_var.set(val_max)
        self.plot_left_histograms()
        self.update_left_display()

    def update_right_red_stretch(self, val_min, val_max):
        self.right_red_min_stretch_var.set(val_min)
        self.right_red_max_stretch_var.set(val_max)
        self.plot_right_histograms()
        self.update_right_display()

    def update_right_green_stretch(self, val_min, val_max):
        self.right_green_min_stretch_var.set(val_min)
        self.right_green_max_stretch_var.set(val_max)
        self.plot_right_histograms()
        self.update_right_display()

    def update_right_blue_stretch(self, val_min, val_max):
        self.right_blue_min_stretch_var.set(val_min)
        self.right_blue_max_stretch_var.set(val_max)
        self.plot_right_histograms()
        self.update_right_display()

    def create_left_histogram_span_selectors(self):
        if not hasattr(self, "left_hist_span_selectors"):
            self.left_hist_span_selectors = []

        if not hasattr(self, "left_hist_figure"):
            return

        for selector in self.left_hist_span_selectors:
            selector.disconnect_events()

        self.left_hist_span_selectors.clear()

        left_red_selector = SpanSelector(
            self.left_hist_axes[0],
            self.update_left_red_stretch,
            "horizontal",
            useblit=True,
        )
        left_green_selector = SpanSelector(
            self.left_hist_axes[1],
            self.update_left_green_stretch,
            "horizontal",
            useblit=True,
        )
        left_blue_selector = SpanSelector(
            self.left_hist_axes[2],
            self.update_left_blue_stretch,
            "horizontal",
            useblit=True,
        )

        self.left_hist_span_selectors.extend(
            [left_red_selector, left_green_selector, left_blue_selector]
        )

    def create_right_histogram_span_selectors(self):
        if not hasattr(self, "right_hist_span_selectors"):
            self.right_hist_span_selectors = []

        if not hasattr(self, "right_hist_figure"):
            return

        for selector in self.right_hist_span_selectors:
            selector.disconnect_events()

        self.right_hist_span_selectors.clear()

        right_red_selector = SpanSelector(
            self.right_hist_axes[0],
            self.update_right_red_stretch,
            "horizontal",
            useblit=True,
        )
        right_green_selector = SpanSelector(
            self.right_hist_axes[1],
            self.update_right_green_stretch,
            "horizontal",
            useblit=True,
        )
        right_blue_selector = SpanSelector(
            self.right_hist_axes[2],
            self.update_right_blue_stretch,
            "horizontal",
            useblit=True,
        )

        self.right_hist_span_selectors.extend(
            [right_red_selector, right_green_selector, right_blue_selector]
        )

    def clear_axes(self):
        # clear left and right axes
        left_xlim = self.left_ax.get_xlim()
        left_ylim = self.left_ax.get_ylim()
        self.left_ax.clear()
        self.left_ax.set_xlim(left_xlim)
        self.left_ax.set_ylim(left_ylim)
        self.left_canvas.draw()

        right_xlim = self.right_ax.get_xlim()
        right_ylim = self.right_ax.get_ylim()
        self.right_ax.clear()
        self.right_ax.set_xlim(right_xlim)
        self.right_ax.set_ylim(right_ylim)
        self.right_canvas.draw()

    # ----------------------------------------------------------------
    # canvas click functionality
    # ----------------------------------------------------------------
    def on_left_canvas_click(self, event):
        # this is called left canvas click but it actually is called for both
        if hasattr(self, "left_data"):
            if (
                self.draw_polygons
                and not self.left_nav_toolbar.mode == "pan/zoom"
                and not self.right_nav_toolbar.mode == "pan/zoom"
            ):
                if event.inaxes == self.left_ax or event.inaxes == self.right_ax:
                    if event.button == 1:
                        # on click and drag, collect all points
                        self.draw_poly_motion_flag = True

                        # # collect one point at a time
                        # x, y = event.xdata, event.ydata
                        # self.current_polygon.append((x, y))
                        # point_color = self.polygon_color_var.get()
                        # self.left_ax.plot(x, y, color = point_color, marker='o')
                        # self.left_canvas.draw()
                        # self.right_ax.plot(x, y, color=point_color, marker='o')
                        # self.right_canvas.draw()
                    elif event.button == 3:
                        if self.template_index is not None:
                            self.best_denom_id_list.append(np.nan)
                            self.ratio_spectrum_list.append(np.nan)
                            if (
                                self.template_index == self.template_index_cache
                                and self.current_denom_polygon_index is not None
                            ):
                                print(self.template_index, self.template_index_cache)
                                # return the item from self.polygon_table at row index self.current_denom_polygon_index
                                # if self.current_denom_polygon_index is None:
                                #     self.current_denom_polygon_index = len(self.all_polygons)
                                item = self.polygon_table_item_id_row[
                                    self.current_denom_polygon_index
                                ]
                                self.polygon_table.delete(item)
                                self.all_polygons.pop(self.current_denom_polygon_index)
                                self.polygon_colors.pop(
                                    self.current_denom_polygon_index
                                )
                                self.best_denom_id_list.pop(
                                    self.current_denom_polygon_index
                                )
                                self.ratio_spectrum_list.pop(
                                    self.current_denom_polygon_index
                                )
                                if self.polygon_spectra:
                                    self.polygon_spectra.pop(
                                        self.current_denom_polygon_index
                                    )
                                if self.polygons_spectral_window is not None:
                                    if self.polygons_spectral_window.winfo_exists():
                                        self.update_polygons_spectral_plot()
                                self.remove_polygons_from_display()
                                self.draw_all_polygons()

                            self.template_index_cache = self.template_index
                            template_polygon = self.all_polygons[self.template_index]
                            denom_x_centroid, denom_y_centroid = (
                                event.xdata,
                                event.ydata,
                            )
                            x_centroid, y_centroid = np.nanmean(
                                template_polygon, axis=0
                            )

                            dx, dy = (
                                denom_x_centroid - x_centroid,
                                denom_y_centroid - y_centroid,
                            )
                            if self.lock_column:
                                print('=== COLUMN LOCK DEBUG ===')
                                print(f'template polygon centroid: ({x_centroid:.1f}, {y_centroid:.1f})')
                                print(f'right-click at: ({denom_x_centroid:.1f}, {denom_y_centroid:.1f})')
                                print(f'initial dx={dx:.1f}, dy={dy:.1f}')
                                # Use the geometry-aware matcher: it picks
                                # the proper backplane (mrrde for tiles,
                                # _in for per-strip), looks up bands by
                                # description, and (when Segment ID is
                                # available) restricts the destination-row
                                # search to the same source strip.
                                locked_dx = self._compute_column_lock_dx(
                                    template_polygon, denom_y_centroid
                                )
                                if locked_dx is not None:
                                    dx = locked_dx
                                else:
                                    # Failure semantics preserved: dx=0
                                    # leaves the polygon at the source x
                                    # rather than jumping to an arbitrary
                                    # click position.
                                    dx = 0
                            print(f'final dx={dx:.1f}, dy={dy:.1f}')
                            self.denom_polygon = [
                                (point[0] + dx, point[1] + dy)
                                for point in template_polygon
                            ]
                            denom_cx, denom_cy = np.nanmean(self.denom_polygon, axis=0)
                            print(f'denom polygon centroid: ({denom_cx:.1f}, {denom_cy:.1f})')
                            print('=== END COLUMN LOCK DEBUG ===')

                            self.clear_axes()
                            self.update_left_display()
                            self.update_right_display()

                            self.current_denom_polygon_index = len(self.all_polygons)
                            self.all_polygons.append(self.denom_polygon)
                            self.polygon_colors.append(self.polygon_color_var.get())
                            self.best_denom_id_list.append(np.nan)
                            self.ratio_spectrum_list.append(np.nan)

                            self.extract_spectra_from_single_polygon()
                            # self.extract_spectra_from_polygons()
                            # self.update_polygons_table()
                            self.draw_all_polygons()

            elif (
                event.inaxes == self.left_ax
                or event.inaxes == self.right_ax
                and not self.left_nav_toolbar.mode == "pan/zoom"
                and not self.right_nav_toolbar.mode == "pan/zoom"
            ):
                if event.button == 1:  # left mouse button press
                    x, y = int(event.xdata), int(event.ydata)
                    self.current_spec_loc = event.xdata, event.ydata
                    self.spectrum_label = f"Pixel: {x}, {y}"
                    if self.bands_inverted:
                        self.spectrum = self.left_data[y, x, :].flatten()
                        self.spectrum = self.spectrum[::-1]
                        self.parameter_values = self.right_data[y, x, :].flatten()
                    else:
                        self.spectrum = self.left_data[y, x, :].flatten()
                        self.parameter_values = self.right_data[y, x, :].flatten()
                    self.spectrum = np.where(self.spectrum < -1, np.nan, self.spectrum)
                    self.spectrum = np.where(self.spectrum > 1.5, np.nan, self.spectrum)
                    self.update_spectral_plot()
                    # clear any previous x's on the image frames
                    self.clear_axes()
                    self.update_left_display()
                    self.update_right_display()

                    # add an x on the image frames where the spectrum was clicked
                    self.left_ax.plot(x, y, "kx")
                    self.right_ax.plot(x, y, "kx")
                    self.left_canvas.draw()
                    self.right_canvas.draw()
                elif event.button == 3:  # Right mouse button press
                    if hasattr(self, "spectrum"):
                        x, y = int(event.xdata), int(event.ydata)

                        # add an x on the image frames where the spectrum was clicked
                        self.left_ax.plot(x, y, "rx")
                        self.right_ax.plot(x, y, "rx")
                        self.left_canvas.draw()
                        self.right_canvas.draw()

                        # plot ratio spectra
                        self.denom_spectrum = self.left_data[y, x, :].flatten()
                        self.denom_spectrum = np.where(
                            self.denom_spectrum < 0, np.nan, self.denom_spectrum
                        )
                        self.denom_spectrum = np.where(
                            self.denom_spectrum > 1.5, np.nan, self.denom_spectrum
                        )
                        self.ratio_spectrum = self.spectrum / self.denom_spectrum
                        self.update_ratio_spectral_plot()
        else:
            messagebox.showwarning(
                "Warning",
                "No data loaded. Load hyperspectral data first into left frame.",
            )

    def on_canvas_motion(self, event):
        if hasattr(self, "left_data") and event.inaxes:
            # Get coordinates in native CRS (typically meters for projected data)
            x_proj, y_proj = self.left_rio.xy(event.ydata, event.xdata)

            # Transform to geographic coordinates (decimal degrees)
            if self.left_rio.crs is not None:
                try:
                    # Determine the appropriate geographic CRS based on the source CRS
                    source_crs = self.left_rio.crs

                    # For Mars data (IAU codes typically 49900+), use IAU:49900
                    # For Earth data (EPSG codes), use EPSG:4326
                    if 'IAU' in str(source_crs).upper() or '49900' in str(source_crs):
                        target_crs = 'IAU:49900'  # Mars geographic
                    else:
                        target_crs = 'EPSG:4326'  # Earth geographic (WGS84)

                    # Transform from projected to geographic coordinates
                    lon, lat = rio_transform(
                        source_crs,
                        target_crs,
                        [x_proj],
                        [y_proj]
                    )
                    coordinate_text = f"Lat: {lat[0]:.6f}, Lon: {lon[0]:.6f}"
                except Exception as e:
                    # If transformation fails, show native coordinates
                    coordinate_text = f"X: {x_proj:.2f}, Y: {y_proj:.2f} (meters)"
            else:
                # Fallback if no CRS info
                coordinate_text = f"X: {x_proj:.2f}, Y: {y_proj:.2f}"

            self.coordinates_label.config(text=coordinate_text)
            self.left_canvas.draw()
        if (
            self.draw_polygons
            and self.draw_poly_motion_flag
            and event.inaxes
            and not self.left_nav_toolbar.mode == "pan/zoom"
        ):
            x, y = event.xdata, event.ydata
            self.current_polygon.append((x, y))
            point_color = self.polygon_color_var.get()
            self.left_ax.plot(x, y, color=point_color, marker=".")
            self.left_canvas.draw()
            self.right_ax.plot(x, y, color=point_color, marker=".")
            self.right_canvas.draw()

    def on_left_release(self, event):
        if hasattr(self, "left_data") and self.left_ax.in_axes(event):
            if self.left_nav_toolbar.mode == "pan/zoom":
                self.left_ax.set_xlim(self.left_ax.get_xlim())
                self.left_ax.set_ylim(self.left_ax.get_ylim())
                self.left_canvas.draw()
                self.right_ax.set_xlim(self.left_ax.get_xlim())
                self.right_ax.set_ylim(self.left_ax.get_ylim())
                self.right_canvas.draw()

                if self.polygons_on_flag:
                    if len(self.all_polygons) > 0:
                        self.draw_all_polygons()

        if self.draw_poly_motion_flag:
            self.draw_poly_motion_flag = False
            if self.current_polygon:
                if len(self.current_polygon) > 2:
                    self.current_polygon.append(
                        self.current_polygon[0]
                    )  # close the polygon
                    polygon_color = self.polygon_color_var.get()

                    self.clear_axes()
                    self.update_left_display()
                    self.update_right_display()

                    self.all_polygons.append(self.current_polygon)
                    self.polygon_colors.append(polygon_color)
                    self.template_index = len(self.all_polygons) - 1
                    self.best_denom_id_list.append(np.nan)
                    self.ratio_spectrum_list.append(np.nan)

                    self.current_polygon = []
                    # self.update_polygons_table()
                    self.extract_spectra_from_single_polygon()
                    # self.extract_spectra_from_polygons()
                    self.draw_all_polygons()
                else:
                    self.current_polygon = []
                    messagebox.showwarning(
                        "Warning", "Polygons must have more than 2 points."
                    )

    def on_right_release(self, event):
        if hasattr(self, "right_data") and self.right_ax.in_axes(event):
            if self.right_nav_toolbar.mode == "pan/zoom":
                self.right_ax.set_xlim(self.right_ax.get_xlim())
                self.right_ax.set_ylim(self.right_ax.get_ylim())
                self.right_canvas.draw()
                self.left_ax.set_xlim(self.right_ax.get_xlim())
                self.left_ax.set_ylim(self.right_ax.get_ylim())
                self.left_canvas.draw()

                if self.polygons_on_flag:
                    if len(self.all_polygons) > 0:
                        self.draw_all_polygons()

        if self.draw_poly_motion_flag:
            self.draw_poly_motion_flag = False
            if self.current_polygon:
                if len(self.current_polygon) > 2:
                    self.current_polygon.append(
                        self.current_polygon[0]
                    )  # close the polygon
                    polygon_color = self.polygon_color_var.get()

                    self.clear_axes()
                    self.update_left_display()
                    self.update_right_display()

                    self.all_polygons.append(self.current_polygon)
                    self.polygon_colors.append(polygon_color)
                    self.template_index = len(self.all_polygons) - 1
                    self.best_denom_id_list.append(np.nan)
                    self.ratio_spectrum_list.append(np.nan)

                    self.current_polygon = []
                    self.extract_spectra_from_single_polygon()
                    # self.extract_spectra_from_polygons()
                    self.draw_all_polygons()
                else:
                    self.current_polygon = []
                    messagebox.showwarning(
                        "Warning", "Polygons must have more than 2 points."
                    )

    def on_scroll(self, event):
        # Handle scroll event for the left canvas
        x_data, y_data = event.xdata, event.ydata
        x_lim, y_lim = self.left_ax.get_xlim(), self.left_ax.get_ylim()

        # Define the zoom factor (adjust as needed)
        zoom_factor = 1.1 if event.button == "down" else 1 / 1.1  # Zoom in or out

        # Adjust the axis limits centered on the mouse cursor position
        new_x_lim = [
            x_data - (x_data - x_lim[0]) * zoom_factor,
            x_data + (x_lim[1] - x_data) * zoom_factor,
        ]
        new_y_lim = [
            y_data - (y_data - y_lim[0]) * zoom_factor,
            y_data + (y_lim[1] - y_data) * zoom_factor,
        ]

        # Update the left canvas
        self.left_ax.set_xlim(new_x_lim)
        self.left_ax.set_ylim(new_y_lim)
        self.left_canvas.draw()

        # Update the right canvas with the same scroll event
        if hasattr(self, "right_canvas"):
            self.right_ax.set_xlim(new_x_lim)
            self.right_ax.set_ylim(new_y_lim)
            self.right_canvas.draw()

    # ----------------------------------------------------------------
    # polygon functions
    # ----------------------------------------------------------------
    def create_polygons_menu_window(self):
        self.polygons_menu_window = tk.Toplevel(self.root)
        self.polygons_menu_window.title("Polygons Menu")

        # ----------------------------------------------------------------
        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.polygons_menu_window)
        ui_frame.grid(row=0, column=0, columnspan=7)

        # ----------------------------------------------------------------
        # create variables to track the row and column position of the buttons
        self.polygons_menu_row = 0
        self.polygons_menu_col = 0

        # ----------------------------------------------------------------
        # Create a button to toggle drawing polygons
        self.draw_polygons_button = tk.Button(
            ui_frame, text="Drawing Polygons Off", command=self.toggle_polygons
        )
        self.draw_polygons_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to toggle column lock
        self.lock_column_button = tk.Button(
            ui_frame, text="Column Lock Off", command=self.toggle_column_lock
        )
        self.lock_column_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to remove polygons from display
        self.remove_polygons_button = tk.Button(
            ui_frame,
            text="Remove Polygons from Display",
            command=self.remove_polygons_from_display,
        )
        self.remove_polygons_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to re-plot the polygons on the display
        self.redraw_all_polygons_button = tk.Button(
            ui_frame, text="Re-draw Polygons on Display", command=self.draw_all_polygons
        )
        self.redraw_all_polygons_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to recalculate parameters for all polygons
        self.recalculate_params_button = tk.Button(
            ui_frame, text="Recalculate Parameters", command=self.recalculate_polygon_parameters
        )
        self.recalculate_params_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to clear all polygons
        self.clear_polygons_button = tk.Button(
            ui_frame, text="Delete All Polygons", command=self.clear_all_polygons
        )
        self.clear_polygons_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to extract spectra from polygons
        self.extract_spectra_button = tk.Button(
            ui_frame,
            text="Plot Mean Spectra",
            command=self.update_polygons_spectral_plot,
        )
        self.extract_spectra_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to save the ROIs
        self.save_rois_button = tk.Button(
            ui_frame, text="Save ROIs", command=self.save_polygons
        )
        self.save_rois_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to load the ROIs
        self.load_rois_button = tk.Button(
            ui_frame, text="Load ROIs", command=self.load_polygons
        )
        self.load_rois_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a dropdown menu to select the polygon color from a list of colors
        self.polygon_color_var = tk.StringVar(ui_frame)
        self.polygon_color_var.set("black")  # default value
        mcolors_keys = mcolors.CSS4_COLORS.keys()
        colors_list = [value for value in mcolors_keys]
        self.polygon_color_dropdown = tk.OptionMenu(
            ui_frame, self.polygon_color_var, *colors_list
        )
        # self.polygon_color_dropdown = tk.OptionMenu(ui_frame, self.polygon_color_var, "red", "green", "blue", "yellow", "orange", "purple", "pink", "brown", "black", "white")
        self.polygon_color_dropdown.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_row += 1

        # ----------------------------------------------------------------
        # plot buttons that act on the currently selected table row (a
        # mouse-friendly alternative to the middle-click context menu)
        self.polygons_menu_col = 0
        self.plot_spectrum_button = tk.Button(
            ui_frame,
            text="Plot Spectrum",
            command=lambda: self.plot_selected_polygon(d=False),
        )
        self.plot_spectrum_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_col += 1

        self.plot_ratio_spectrum_button = tk.Button(
            ui_frame,
            text="Plot Ratio Spectrum",
            command=lambda: self.plot_selected_polygon(d=True),
        )
        self.plot_ratio_spectrum_button.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_row += 1

        # ----------------------------------------------------------------
        # create a table to display polygon information
        self.create_polygons_table()

    # ------------------------------------------------
    # polygon table functions
    # ------------------------------------------------
    def create_polygons_table(self):
        self.polygon_table_header_labels = (
            "Polygon Number",
            "Color",
            "Number of Points",
            "Denominator",
            "Template",
            "Mineral ID 1",
            "Mineral ID 2",
            "Mineral ID 3",
            "Mineral ID 4",
            "wvl",
            "Spectrum Mean",
            "params",
            "Parameters Mean",
            "Best Denom ID",
            "Ratio Spectrum",
        )
        self.polygon_editing_data = {}
        self.editing_polygon_table = False

        # Create a Text widget for displaying the header labels
        self.polygon_table = ttk.Treeview(
            self.polygons_menu_window,
            columns=self.polygon_table_header_labels,
            show="headings",
        )

        # Configure column headings
        for col in self.polygon_table_header_labels:
            self.polygon_table.heading(col, text=col)
            self.polygon_table.column(
                col, width=len(col) * 9, anchor=tk.CENTER
            )  # Set the column width as needed

        self.polygon_table.grid(
            row=1, column=0, columnspan=len(self.polygon_table_header_labels)
        )

        # Bind right-click to show the context menu
        self.polygon_table.bind("<Button-2>", self.show_context_menu)

        # bind single-click to check a box
        self.polygon_table.bind("<Button-1>", self.check_a_box)

        # highlight selected polygon
        self.polygon_table.bind("<ButtonRelease-1>", self.highlight_polygon)

        # double click to edit a cell
        self.polygon_table.bind("<Double-1>", self.edit_cell)

        # Create a context menu
        self.context_menu = tk.Menu(self.polygon_table, tearoff=0)
        self.context_menu.add_command(label="Edit", command=self.edit_cell)
        self.context_menu.add_command(label="Delete", command=self.delete_polygon)
        self.context_menu.add_command(
            label="Plot Spectrum", command=self.plot_selected_polygon
        )
        self.context_menu.add_command(
            label="Plot Ratio Spectrum",
            command=lambda: self.plot_selected_polygon(d=True),
        )

        self.update_polygons_table()

    def update_polygons_table(self):
        # get the information from the polygons to display in the table 13 14
        if hasattr(self, "polygons_table_data"):
            if len(self.polygons_table_data) != 0:
                min_ids = [line[5:9] for line in self.polygons_table_data]
                if self.upload_polygons_flag == False:
                    try:
                        self.best_denom_id_list = [
                            line[13] for line in self.polygons_table_data
                        ]
                        self.ratio_spectrum_list = [
                            line[14] for line in self.polygons_table_data
                        ]
                    except:
                        print("No best denom id or ratio spectrum")
        self.polygons_table_data = []
        for i, polygon in enumerate(self.all_polygons):
            polygon_number = i
            polygon_color = self.polygon_colors[i]
            number_of_points = len(polygon) - 1
            denom_check = "[ ]"
            template_check = "[ ]"
            if self.add_spectrum_flag and self.num_idx == i:
                if self.denominator_index is not None:
                    best_denom_id = self.denominator_index
                    row_ratio_spectrum = self.spectrum
            else:
                try:
                    best_denom_id = self.best_denom_id_list[i]
                except:
                    best_denom_id = np.nan
                try:
                    row_ratio_spectrum = self.ratio_spectrum_list[i]
                except:
                    row_ratio_spectrum = np.nan
            if self.denominator_index == i:
                denom_check = "[X]"
            if self.template_index == i:
                template_check = "[X]"
            if "min_ids" in locals() and i < len(min_ids):
                minid1, minid2, minid3, minid4 = min_ids[i]
            else:
                minid1, minid2, minid3, minid4 = "", "", "", ""
            try:
                spec_mean = self.polygon_spectra[i]
                param_mean = self.polygon_params[i]
            except IndexError:
                print(f"No means for polygon {i}")

            row = (
                polygon_number,
                polygon_color,
                number_of_points,
                denom_check,
                template_check,
                minid1,
                minid2,
                minid3,
                minid4,
                self.left_wvl,
                spec_mean,
                self.right_wvl,
                param_mean,
                best_denom_id,
                row_ratio_spectrum,
            )
            self.polygons_table_data.append(row)

        for item in self.polygon_table.get_children():
            self.polygon_table.delete(item)

        self.polygon_table_item_id_row = []
        for row in self.polygons_table_data:
            item_id = self.polygon_table.insert("", "end", values=row)
            self.polygon_table_item_id_row.append(item_id)
        self.add_spectrum_flag = False
        self.upload_polygons_flag = False

    def highlight_polygon(self, event):
        item = self.polygon_table.identify_row(event.y)
        if item:
            self.polygon_table.selection_set(item)
            # get the polygon number
            self.polygon_to_highlight = int(self.polygon_table.item(item, "values")[0])
            self.draw_all_polygons()

    def check_a_box(self, event):
        if event is None:
            event = self.tmp_event
        item = self.polygon_table.identify_row(event.y)
        column = self.polygon_table.identify_column(event.x)
        if int(column[1:]) != 4 and int(column[1:]) != 5:
            return
        cell_value = self.polygon_table.item(item, "values")[int(column[1:]) - 1]

        if int(column[1:]) == 4:
            if cell_value == "[ ]":
                # set all the other items in the column to '[ ]'
                for item_ in self.polygon_table.get_children():
                    self.polygon_table.set(item_, column, "[ ]")
                self.polygon_table.set(item, column, "[X]")
                # populate self.denominator index with the value of polygon Number column
                self.denominator_index = int(self.polygon_table.item(item, "values")[0])
            else:
                self.polygon_table.set(item, column, "[ ]")
                self.denominator_index = None

        if int(column[1:]) == 5:
            if cell_value == "[ ]":
                # set all the other items in the column to '[ ]'
                for item_ in self.polygon_table.get_children():
                    self.polygon_table.set(item_, column, "[ ]")
                self.polygon_table.set(item, column, "[X]")
                # populate self.denominator index with the value of polygon Number column
                self.template_index = int(self.polygon_table.item(item, "values")[0])
            else:
                self.polygon_table.set(item, column, "[ ]")
                self.template_index = []
                self.template_index_cache = []
            # Source strip changed → refresh the footprint outline.
            self._update_strip_footprint_overlay()

    def clear_polygons_table(self):
        self.polygons_table_data = []
        for i, polygon in enumerate(self.all_polygons):
            polygon_number = None
            polygon_color = None
            number_of_points = None
            self.polygons_table_data.append(
                [polygon_number, polygon_color, number_of_points]
            )

        # Populate the table with data
        for i, row_data in enumerate(self.polygons_table_data):
            for j, cell_value in enumerate(row_data):
                cell_label = tk.Label(self.polygons_menu_window, text=cell_value)
                cell_label.grid(row=i + 2, column=j)

    def show_context_menu(self, event):
        self.tmp_event = event
        item = self.polygon_table.identify_row(event.y)
        if item:
            self.polygon_table.selection_set(item)
            self.context_menu.post(event.x_root, event.y_root)

    def edit_cell(self, event=None):
        # Store event for later use (right-click menu calls with event=None)
        if event is not None:
            self.tmp_event = event
        event = self.tmp_event

        # Get item and column from event coordinates
        item = self.polygon_table.identify_row(event.y)
        column = self.polygon_table.identify_column(event.x)

        # Fallback to current selection if identify_row returns empty
        if not item:
            selection = self.polygon_table.selection()
            if selection:
                item = selection[0]
            else:
                return

        # Check if column is valid and editable
        if not column or column == '#0':
            return
        col_num = int(column[1:])
        if col_num not in (2, 6, 7, 8, 9):
            return

        cell_value = self.polygon_table.item(item, "values")[col_num - 1]
        # Get polygon number (first column) - this is stable across table refreshes
        pc_index = int(self.polygon_table.item(item, "values")[0])

        if not self.editing_polygon_table:
            self.polygon_editing_data[item, column] = (
                cell_value  # Store the original value
            )

            # Create a Toplevel window for editing
            edit_window = tk.Toplevel(self.polygons_menu_window)

            # Create an Entry widget inside the edit_window
            edit_entry = ttk.Entry(edit_window, justify="center")
            edit_entry.insert(0, cell_value)
            edit_entry.pack()

            # Capture values for closure - use pc_index instead of item ID
            # since item IDs can change if the table is refreshed
            captured_pc_index = pc_index
            captured_column = column
            captured_col_num = col_num

            def save_and_close():
                new_value = edit_entry.get()

                # Find the current item ID for this polygon number
                current_item = None
                for tree_item in self.polygon_table.get_children():
                    if int(self.polygon_table.item(tree_item, "values")[0]) == captured_pc_index:
                        current_item = tree_item
                        break

                if current_item is None:
                    print(f"Warning: Could not find polygon {captured_pc_index} in table")
                    edit_window.destroy()
                    return

                self.polygon_table.set(current_item, captured_column, new_value)
                if captured_col_num == 2:
                    self.polygon_colors[captured_pc_index] = new_value
                else:
                    updated_row = list(self.polygons_table_data[captured_pc_index])
                    updated_row[captured_col_num - 1] = new_value
                    self.polygons_table_data[captured_pc_index] = tuple(updated_row)
                self.draw_all_polygons()
                self.update_polygons_spectral_plot()
                edit_window.destroy()

            # Create a button to save changes
            save_button = ttk.Button(edit_window, text="Save", command=save_and_close)
            save_button.pack()

            # Bind the "Return" key to the "Save" button's functionality
            edit_window.bind("<Return>", lambda e: save_and_close())

            # Focus on the Entry widget
            edit_entry.focus_set()

            # self.update_polygons_spectral_plot()

    def add_spectrum_to_table(self):
        self.add_spectrum_flag = True
        self.update_polygons_table()

    # ------------------------------------------------
    # polygon menu functions
    # ------------------------------------------------
    def toggle_polygons(self):
        if self.draw_polygons:
            self.draw_polygons = False
            self.current_polygon = []
            self.draw_polygons_button.config(
                text="Drawing Polygons Off", relief="sunken"
            )
        else:
            self.draw_polygons = True
            self.draw_polygons_button.config(
                text="Drawing Polygons On", relief="raised"
            )

    def toggle_column_lock(self):
        if self.lock_column:
            self.lock_column = False
            self.lock_column_button.config(text="Lock Column Off", relief="sunken")
        else:
            self.lock_column = True
            self.lock_column_button.config(text="Lock Column On", relief="raised")
        # Show/hide the source-strip footprint outline.
        self._update_strip_footprint_overlay()

    def _get_column_lock_geometry(self):
        """Return cached IR Sample / Segment ID arrays for the current cube.

        On first call for a given source file, resolves the proper
        geometry companion (mrrde for tiles, _in for per-strip
        products), opens it with rasterio, and caches the IR Sample
        and (when present) Segment ID bands looked up *by description*
        rather than hardcoded index — band 6 = IR Sample is only true
        for per-strip _in files; tile mrrde files put IR Sample at
        band 5.

        Returns a dict with keys `ir_sample`, `segment_id` (may be
        None on per-strip products that don't carry it), `nodata`,
        and `path`. Returns None if no geometry file is available —
        the caller must then disable column-lock for this cube and
        tell the user, instead of falling back to reading the source
        as a fake backplane (which is what produced the random-paste
        bug on tile data).

        Negative results (no geometry file, unreadable file, missing
        IR Sample band) are cached as None so the regex/file probe
        doesn't re-run on every paste.
        """
        source_file = self.left_data.filename
        if source_file in self._geometry_cache:
            return self._geometry_cache[source_file]

        geom_path = _resolve_geometry_file(source_file)
        if geom_path is None:
            self._geometry_cache[source_file] = None
            return None

        try:
            with rio.open(geom_path) as geom:
                ir_idx = _find_band_index(geom, 'ir_sample')
                tgt_idx = _find_band_index(geom, 'target_id')
                seg_idx = _find_band_index(geom, 'segment_id')
                if ir_idx is None:
                    print(
                        f'column-lock: geometry file {geom_path} has no '
                        f'IR Sample band; descriptions: {geom.descriptions}'
                    )
                    self._geometry_cache[source_file] = None
                    return None
                ir_sample = geom.read(ir_idx).astype(np.float64)
                target_id = geom.read(tgt_idx).astype(np.float64) if tgt_idx else None
                segment_id = geom.read(seg_idx).astype(np.float64) if seg_idx else None
                nodata = geom.nodata
        except Exception as e:
            print(f'column-lock: failed to read geometry file {geom_path}: {e}')
            self._geometry_cache[source_file] = None
            return None

        # Treat the file's nodata sentinel and NaN as missing across all
        # arrays; downstream code can then rely on np.isnan() uniformly.
        if nodata is not None:
            ir_sample = np.where(ir_sample == nodata, np.nan, ir_sample)
            if target_id is not None:
                target_id = np.where(target_id == nodata, np.nan, target_id)
            if segment_id is not None:
                segment_id = np.where(segment_id == nodata, np.nan, segment_id)

        entry = {
            'ir_sample': ir_sample,
            'target_id': target_id,
            'segment_id': segment_id,
            'nodata': nodata,
            'path': geom_path,
        }
        self._geometry_cache[source_file] = entry
        return entry

    def _compute_column_lock_dx(self, template_polygon, click_y):
        """Compute dx for a column-locked paste.

        The template polygon defines a sensor column (mean IR Sample
        under the polygon's filled mask) and, when the geometry file
        carries them, a *source strip identity*:

          - Target ID (band 1) — the unique CRISM observation ID;
            this is the actual strip-identifier in mosaic tiles.
          - Segment ID (band 2) — a within-target counter; only useful
            as a secondary gate for multi-segment observations.

        At the destination row given by `click_y`, the search is
        gated to pixels that share the source's Target ID and (if
        present) Segment ID, then picks the one with IR Sample
        closest to the source mean.

        Why both: an earlier version gated only on Segment ID and
        still produced wrong-strip placements. Segment ID takes only
        a handful of distinct values across a mosaic (e.g. {1,3,5,7}),
        so many physically distinct strips share the same value.
        Target ID is the strip-unique key.

        Returns the dx (relative to the template's centroid x) that
        places the polygon's centroid at the matched x. Returns None
        when column-lock cannot be applied for this paste — most
        commonly because the source strip simply doesn't reach the
        destination row the user clicked. Hard failures (no geometry
        file) raise a one-shot messagebox per source cube.
        """
        geom = self._get_column_lock_geometry()
        if geom is None:
            source_file = self.left_data.filename
            if source_file not in self._geometry_warning_shown:
                self._geometry_warning_shown.add(source_file)
                messagebox.showwarning(
                    "Column lock",
                    "No CRISM geometry/backplane file found for "
                    f"{os.path.basename(source_file)}.\n\n"
                    "Column-lock cannot determine sensor columns "
                    "without a backplane file (e.g. _mrrde_ for "
                    "tiles, _in*_ for per-strip products) and will "
                    "be skipped for this cube."
                )
            return None

        ir_sample = geom['ir_sample']
        target_id = geom['target_id']
        segment_id = geom['segment_id']
        nrows, ncols = ir_sample.shape

        # Compute source identity over the polygon's filled footprint
        # (more representative than just its vertices for large ROIs).
        mask = np.zeros((nrows, ncols), dtype=np.uint8)
        pts = np.array([[int(x), int(y)] for x, y in template_polygon])
        cv2.fillPoly(mask, [pts], 1)
        mask_bool = mask.astype(bool)

        src_ir_vals = ir_sample[mask_bool]
        src_ir_vals = src_ir_vals[~np.isnan(src_ir_vals)]
        if src_ir_vals.size == 0:
            print('column-lock: template polygon has no valid IR Sample values')
            return None
        src_ir_mean = float(np.mean(src_ir_vals))

        def _mode_under_mask(arr):
            if arr is None:
                return None
            vals = arr[mask_bool]
            vals = vals[~np.isnan(vals)]
            if vals.size == 0:
                return None
            # Mode, not mean: averaging strip identifiers makes no
            # physical sense; a polygon straddling two strips picks
            # whichever covers more of the polygon.
            uniq, counts = np.unique(vals.astype(np.int64), return_counts=True)
            return int(uniq[int(np.argmax(counts))])

        src_target = _mode_under_mask(target_id)
        src_seg = _mode_under_mask(segment_id)

        dest_y = int(click_y)
        if not (0 <= dest_y < nrows):
            print(f'column-lock: destination row y={dest_y} out of bounds (0..{nrows-1})')
            return None

        row_ir = ir_sample[dest_y, :]
        valid = ~np.isnan(row_ir)
        if src_target is not None:
            valid &= (target_id[dest_y, :] == src_target)
        if src_seg is not None and segment_id is not None:
            valid &= (segment_id[dest_y, :] == src_seg)
        valid_xs = np.where(valid)[0]

        if valid_xs.size == 0:
            # The source strip simply doesn't reach this destination
            # row. This is the most common reason column-lock can't
            # apply: users click somewhere along-track that's outside
            # the source observation's footprint within the mosaic.
            if src_target is not None:
                print(
                    f'column-lock: source strip (target {src_target}'
                    + (f', segment {src_seg}' if src_seg is not None else '')
                    + f') does not reach destination row y={dest_y}; '
                    'cannot lock — falling back to no horizontal shift'
                )
            else:
                print(f'column-lock: no valid IR Sample pixels in destination row y={dest_y}')
            return None

        diffs = np.abs(row_ir[valid_xs] - src_ir_mean)
        best_x = int(valid_xs[int(np.argmin(diffs))])

        # Source polygon's centroid x (independent of the click).
        src_x = float(np.nanmean([pt[0] for pt in template_polygon]))

        print(
            f'column-lock: src target={src_target}, segment={src_seg}, '
            f'IR Sample={src_ir_mean:.2f} → '
            f'dest x={best_x}, IR Sample={row_ir[best_x]:.2f}, '
            f'dx={best_x - src_x:.1f}'
        )
        return best_x - src_x

    def _source_strip_identity(self):
        """Return (target_id, segment_id) for the current template polygon.

        Uses the *mode* of the geometry bands under the template's
        filled footprint, so a polygon that straddles a strip boundary
        picks whichever strip covers more of it. Returns (None, None)
        if column-lock can't be applied at all (no template, no
        geometry, no IR Sample band, etc.) — caller treats that as
        "nothing to highlight".
        """
        # `template_index` can be None, [], or a non-negative int
        # (existing convention in this codebase). Treat anything that
        # isn't a usable integer index as "no template selected".
        tpl = getattr(self, 'template_index', None)
        if tpl is None or isinstance(tpl, (list, tuple)):
            return (None, None)
        try:
            tpl_idx = int(tpl)
        except (TypeError, ValueError):
            return (None, None)
        if tpl_idx < 0 or tpl_idx >= len(self.all_polygons):
            return (None, None)

        geom = self._get_column_lock_geometry()
        if geom is None:
            return (None, None)

        target_id = geom['target_id']
        segment_id = geom['segment_id']
        if target_id is None and segment_id is None:
            return (None, None)

        nrows, ncols = geom['ir_sample'].shape
        mask = np.zeros((nrows, ncols), dtype=np.uint8)
        pts = np.array([[int(x), int(y)] for x, y in self.all_polygons[tpl_idx]])
        cv2.fillPoly(mask, [pts], 1)
        mb = mask.astype(bool)

        def _mode(arr):
            if arr is None:
                return None
            v = arr[mb]
            v = v[~np.isnan(v)]
            if v.size == 0:
                return None
            u, c = np.unique(v.astype(np.int64), return_counts=True)
            return int(u[int(np.argmax(c))])

        return (_mode(target_id), _mode(segment_id))

    def _compute_strip_footprint_contours(self, src_target, src_segment):
        """Return list of polygon-vertex arrays outlining the source strip.

        `src_target` and `src_segment` may each be None if the
        corresponding band wasn't present in the geometry file. The
        mask is the AND of whichever identifiers are available; if
        both are None, returns an empty list (no useful outline).
        Cached by (target, segment) to avoid recomputing every redraw.
        """
        key = (src_target, src_segment)
        if key in self._strip_footprint_cache:
            return self._strip_footprint_cache[key]

        geom = self._get_column_lock_geometry()
        if geom is None:
            return []

        target_id = geom['target_id']
        segment_id = geom['segment_id']
        if src_target is None and src_segment is None:
            return []

        # Build the strip mask. A pixel belongs to the source strip
        # iff every available identifier matches.
        if target_id is not None and src_target is not None:
            mask = (target_id == src_target)
            if segment_id is not None and src_segment is not None:
                mask &= (segment_id == src_segment)
        elif segment_id is not None and src_segment is not None:
            mask = (segment_id == src_segment)
        else:
            return []

        mask_u8 = mask.astype(np.uint8)
        # CHAIN_APPROX_TC89_L1 collapses straight-line vertices, which
        # keeps the artist count low without losing visible detail at
        # the typical zoom levels users work at.
        contours, _ = cv2.findContours(
            mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_TC89_L1
        )
        # cv2 returns Nx1x2 int32 arrays of (x, y); flatten to Nx2.
        polys = [c.reshape(-1, 2) for c in contours if len(c) >= 3]
        self._strip_footprint_cache[key] = polys
        return polys

    def _clear_strip_footprint_overlay(self):
        """Remove any drawn footprint outlines from both canvases."""
        for artist in self._strip_footprint_artists:
            try:
                artist.remove()
            except (ValueError, NotImplementedError):
                # Artist may already be gone if its axes were cleared.
                pass
        self._strip_footprint_artists = []

    def _update_strip_footprint_overlay(self, redraw=True):
        """Refresh the source-strip outline on left and right canvases.

        Drawn only when column-lock is on AND a template polygon is
        selected AND the cube has usable geometry. Otherwise removes
        any previous outline. Idempotent — safe to call from any
        place that might change the overlay's visibility (toggle,
        template change, redraw, cube load).
        """
        self._clear_strip_footprint_overlay()
        if not getattr(self, 'lock_column', False):
            if redraw:
                self.left_canvas.draw_idle()
                self.right_canvas.draw_idle()
            return
        if not hasattr(self, 'left_data'):
            return

        src_target, src_segment = self._source_strip_identity()
        if src_target is None and src_segment is None:
            if redraw:
                self.left_canvas.draw_idle()
                self.right_canvas.draw_idle()
            return

        polys = self._compute_strip_footprint_contours(src_target, src_segment)
        if not polys:
            if redraw:
                self.left_canvas.draw_idle()
                self.right_canvas.draw_idle()
            return

        # Yellow outline contrasts with both the polygon ink (black/
        # cyan) and CRISM imagery; double-line (background black,
        # foreground yellow) keeps it readable on light terrain too.
        for ax in (self.left_ax, self.right_ax):
            for verts in polys:
                # Background stroke (slightly wider, dark) for legibility.
                bg = Polygon(
                    verts,
                    closed=True,
                    fill=False,
                    edgecolor='black',
                    linewidth=2.4,
                    zorder=4,
                )
                fg = Polygon(
                    verts,
                    closed=True,
                    fill=False,
                    edgecolor='yellow',
                    linewidth=1.2,
                    zorder=5,
                )
                ax.add_patch(bg)
                ax.add_patch(fg)
                self._strip_footprint_artists.extend([bg, fg])

        if redraw:
            self.left_canvas.draw_idle()
            self.right_canvas.draw_idle()

    def draw_all_polygons(self):
        self.polygons_on_flag = True
        for polygon_color, polygon in zip(self.polygon_colors, self.all_polygons):
            # x, y = zip(*polygon)
            # self.left_ax.plot(x, y, 'ro')
            # self.right_ax.plot(x, y, 'ro')
            if self.highlight_polygon:
                if self.all_polygons.index(polygon) == self.polygon_to_highlight:
                    self.left_ax.add_patch(
                        Polygon(
                            polygon,
                            closed=True,
                            facecolor=polygon_color,
                            edgecolor="cyan",
                        )
                    )
                    self.right_ax.add_patch(
                        Polygon(
                            polygon,
                            closed=True,
                            facecolor=polygon_color,
                            edgecolor="cyan",
                        )
                    )
                else:
                    self.left_ax.add_patch(
                        Polygon(
                            polygon, closed=True, facecolor=polygon_color, edgecolor="k"
                        )
                    )
                    self.right_ax.add_patch(
                        Polygon(
                            polygon, closed=True, facecolor=polygon_color, edgecolor="k"
                        )
                    )
            else:
                self.left_ax.add_patch(
                    Polygon(
                        polygon, closed=True, facecolor=polygon_color, edgecolor="k"
                    )
                )
                self.right_ax.add_patch(
                    Polygon(
                        polygon, closed=True, facecolor=polygon_color, edgecolor="k"
                    )
                )
        # Re-attach the source-strip footprint outline so it survives
        # the polygon redraw. Use redraw=False to fold it into the
        # single canvas.draw() below instead of triggering another.
        self._update_strip_footprint_overlay(redraw=False)
        self.left_canvas.draw()
        self.right_canvas.draw()
        self.update_polygons_table()

    def remove_polygons_from_display(self):
        self.polygons_on_flag = False
        self.clear_axes()
        self.update_left_display()
        self.update_right_display()

    def recalculate_polygon_parameters(self):
        self.polygon_spectra = []
        self.polygon_params = []
        if hasattr(self, "left_data"):
            self.extract_spectra_from_polygons()
            self.update_polygons_spectral_plot()


    def delete_polygon(self):
        event = self.tmp_event
        # Ask the user if they are sure
        answer = messagebox.askyesno(
            "Confirmation", "Are you sure you want to delete this polygon?"
        )
        if answer:
            item = self.polygon_table.identify_row(event.y)
            pc_index = int(self.polygon_table.item(item, "values")[0])
            if len(self.polygons_table_data) != 0:
                self.polygons_table_data.pop(pc_index)
            self.polygon_table.delete(item)
            self.all_polygons.pop(pc_index)
            self.polygon_colors.pop(pc_index)
            self.best_denom_id_list.pop(pc_index)
            self.ratio_spectrum_list.pop(pc_index)
            if self.polygon_spectra:
                self.polygon_spectra.pop(pc_index)
            if self.polygons_spectral_window is not None:
                if self.polygons_spectral_window.winfo_exists():
                    self.update_polygons_spectral_plot()
            self.remove_polygons_from_display()
            self.draw_all_polygons()
            self.current_denom_polygon_index = None

        else:
            pass

    def _selected_polygon_index(self):
        """Return the polygon number of the active table row, or None.

        Prefers self.polygon_to_highlight, the app's persistent "active
        polygon" set when the user clicks a table row (which also highlights
        the polygon on the image). The Treeview's own .selection() is
        unreliable here because update_polygons_table() rebuilds every row,
        which clears it; polygon_to_highlight survives that rebuild. Falls
        back to the live selection, then to the last context-menu click.
        """
        candidates = []

        # Authoritative: set by highlight_polygon on row click. Initialized to
        # False, and polygon number 0 is valid, so test identity, not truth.
        pth = self.polygon_to_highlight
        if pth is not False and pth is not None:
            candidates.append(int(pth))

        # Live Treeview selection (set on the same click, before any rebuild).
        selection = self.polygon_table.selection()
        if selection:
            candidates.append(int(self.polygon_table.item(selection[0], "values")[0]))

        # Last resort: the row under the most recent context-menu click.
        if getattr(self, "tmp_event", None) is not None:
            item = self.polygon_table.identify_row(self.tmp_event.y)
            if item:
                candidates.append(int(self.polygon_table.item(item, "values")[0]))

        # Drop stale indices (e.g. a highlight left over after a deletion).
        for idx in candidates:
            if 0 <= idx < len(self.all_polygons):
                return idx
        return None

    def plot_selected_polygon(self, d=None):
        pc_index = self._selected_polygon_index()
        if pc_index is None:
            messagebox.showwarning(
                "No polygon selected",
                "Click a polygon row in the table first (it will highlight on "
                "the image), then click the plot button.",
            )
            return

        ssw = spectral_window  # single_spectrum_window
        ssw.create_spectral_plot(app, pc_index, d=d)

    def clear_all_polygons(self):
        # Ask the user if they are sure
        answer = messagebox.askyesno(
            "Confirmation", "Are you sure you want to delete all polygons?"
        )
        if answer:
            self.clear_polygons_table()
            self.all_polygons = []
            self.polygon_colors = []
            self.polygon_spectra = []
            self.update_polygons_spectral_plot()
            self.remove_polygons_from_display()
            self.update_polygons_table()
        else:
            pass

    def extract_spectra_from_single_polygon(self, polygon_index=-1):
        if hasattr(self, "left_data"):
            nrows, ncols = self.left_data.nrows, self.left_data.ncols

            # Prepare mask for the specified polygon
            polygon = self.all_polygons[polygon_index]
            mask = np.zeros((nrows, ncols), dtype=np.uint8)
            polygon_points_int = [(int(x), int(y)) for x, y in polygon]
            cv2.fillPoly(mask, [np.array(polygon_points_int)], 1)

            # ROI mean. Filter sentinel/out-of-range *per pixel* before
            # averaging — otherwise a CRISM tile's 65535 fill values
            # (data-ignore sentinel) inside a partial-coverage ROI drag
            # the band means above 1.5 and the post-mean filter wipes
            # the whole spectrum to NaN.
            left_fill = _data_ignore_value(self.left_data)
            right_fill = _data_ignore_value(self.right_data) if hasattr(self, 'right_data') else None
            spec_fills = tuple(v for v in (left_fill, 65535.0) if v is not None)
            param_fills = tuple(v for v in (right_fill, 65535.0) if v is not None)

            def calculate_stats(data, mask, fills, valid_range):
                return _nanmean_over_mask(
                    data, mask, fill_values=fills, valid_range=valid_range
                )

            mean_spec = calculate_stats(self.ld, mask, spec_fills, (0.0, 1.5))
            mean_param = calculate_stats(self.rd, mask, param_fills, None)

            # Calculate statistics for the left data
            if polygon_index != -1:
                print("here")
                self.polygon_spectra[polygon_index] = mean_spec
                self.polygon_params[polygon_index] = mean_param
            else:
                self.polygon_spectra.append(mean_spec)
                self.polygon_params.append(mean_param)

            # Plot the spectra
            self.update_polygons_spectral_plot()
        else:
            messagebox.showwarning(
                "Warning",
                "No data loaded. Load hyperspectral data first into left frame.",
            )

    def extract_spectra_from_polygons(self):
        if hasattr(self, "left_data"):
            self.polygon_spectra = []
            self.polygon_params = []

            if self.all_polygons:
                nrows, ncols = self.left_data.nrows, self.left_data.ncols

                # Prepare masks for all polygons
                print(f"Extracting means for {len(self.all_polygons)} polygons")
                masks = np.zeros((len(self.all_polygons), nrows, ncols), dtype=np.uint8)
                for idx, polygon in enumerate(self.all_polygons):
                    polygon_points_int = [(int(x), int(y)) for x, y in polygon]
                    cv2.fillPoly(masks[idx], [np.array(polygon_points_int)], 1)

                # See note on the single-polygon variant: pre-filter
                # sentinel/out-of-range pixels per band before averaging
                # so partial-coverage ROIs over CRISM tile no-data
                # regions still produce a real spectrum.
                left_fill = _data_ignore_value(self.left_data)
                right_fill = _data_ignore_value(self.right_data) if hasattr(self, 'right_data') else None
                spec_fills = tuple(v for v in (left_fill, 65535.0) if v is not None)
                param_fills = tuple(v for v in (right_fill, 65535.0) if v is not None)

                def calculate_stats(data, mask, fills, valid_range):
                    return _nanmean_over_mask(
                        data, mask, fill_values=fills, valid_range=valid_range
                    )

                # Calculate statistics in parallel
                results = Parallel(n_jobs=-1)(
                    delayed(calculate_stats)(self.ld, mask, spec_fills, (0.0, 1.5))
                    for mask in masks
                )
                print(len(results))
                self.polygon_spectra.extend(results)

                # Calculate parameters (assuming self.right_data is defined)
                results_params = Parallel(n_jobs=-1)(
                    delayed(calculate_stats)(self.rd, mask, param_fills, None)
                    for mask in masks
                )
                self.polygon_params.extend(results_params)

                # plot the spectra
                self.update_polygons_spectral_plot()
            else:
                messagebox.showwarning(
                    "Warning", "No polygons drawn. Draw polygons first."
                )
        else:
            messagebox.showwarning(
                "Warning",
                "No data loaded. Load hyperspectral data first into left frame.",
            )

        # old, slow method
        # if hasattr(self, "left_data"):
        #     if self.all_polygons:

        #         self.polygon_spectra = []
        #         for polygon in self.all_polygons:
        #             mask = np.zeros((self.left_data.nrows, self.left_data.ncols), dtype=np.uint8)
        #             polygon_points_int = [(int(x), int(y)) for x, y in polygon]
        #             cv2.fillPoly(mask, [np.array(polygon_points_int)], 1)
        #             gstats = spectral.calc_stats(self.left_data, mask=mask, allow_nan=False)
        #             pstats = spectral.calc_stats(self.right_data, mask=mask, allow_nan=True)
        #             mean_spectrum = gstats.mean
        #             # set values <=0 or >=1 to np.nan
        #             mean_spectrum = np.where(mean_spectrum < 0, np.nan, mean_spectrum)
        #             mean_spectrum = np.where(mean_spectrum > 1, np.nan, mean_spectrum)

        #             mean_param = pstats.mean

        #             self.polygon_spectra.append(mean_spectrum)
        #             self.polygon_params.append(mean_param)

        #         # plot the spectra
        #         self.update_polygons_spectral_plot()
        #     else:
        #         messagebox.showwarning("Warning", "No polygons drawn. Draw polygons first.")
        # else:
        #     messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first into left frame.")

    def save_polygons(self):
        if self.all_polygons:
            # get the filename to save the polygons to
            filename = filedialog.asksaveasfilename(
                initialdir="/", title="Select file"
            )  # ,filetypes = (("shp files","*.shp")))
            if filename:
                # Add .gpkg extension if not already present
                if not filename.lower().endswith(".gpkg"):
                    filename += ".gpkg"

                # if the file is already present, ask the user if they want to save
                # Check if file already exists
                if os.path.exists(filename):
                    # Ask user if they want to overwrite the file
                    overwrite = messagebox.askyesno(
                        "Overwrite Confirmation",
                        f"The file '{filename}' already exists. Do you want to overwrite it?",
                    )
                    if not overwrite:
                        return
                # save the polygons using geopandas and shapely.geometry
                all_geo_poly_coords = []
                for poly_coords in self.all_polygons:
                    geo_poly_coords = [self.left_rio.xy(y, x) for x, y in poly_coords]
                    all_geo_poly_coords.append(geo_poly_coords)
                polygon_geometries = [
                    sgp(poly_coords) for poly_coords in all_geo_poly_coords
                ]

                # create a GeoDataFrame
                column_dict = {}
                # Assuming each column has a unique identifier (replace 'column_id' with the actual identifier)
                for i, column_id in enumerate(
                    self.polygon_table_header_labels
                ):  # Add the actual identifiers for your columns
                    column_values = [
                        self.polygon_table.item(item_id, "values")[i]
                        for item_id in self.polygon_table.get_children()
                    ]
                    column_dict[column_id] = column_values
                print(column_dict)
                polygons_gdf = gpd.GeoDataFrame(
                    column_dict, geometry=polygon_geometries, crs=self.left_rio.crs
                )

                # polygons_gdf = gpd.GeoDataFrame({'color': self.polygon_colors}, geometry=polygon_geometries, crs=self.left_rio.crs)
                polygons_gdf.to_file(filename, driver="GPKG")
        else:
            messagebox.showwarning("Warning", "No polygons drawn. Draw polygons first.")

    def load_polygons(self):
        self.upload_polygons_flag = True
        self.best_denom_id_list = []
        self.ratio_spectrum_list = []

        # Function to parse and clean a parameter values vector string
        def parse_parameter_values(param_str):
            param_str = param_str.replace("[", "").replace("]", "").replace("\n", "")
            parameter_values = param_str.split()
            parameter_values = np.array([float(val) for val in parameter_values])
            return parameter_values

        # Function to parse and clean a spectrum string
        def parse_spectrum(spectrum_str):
            spectrum_str = (
                spectrum_str.replace("[", "").replace("]", "").replace("\n", "")
            )
            spectrum_values = spectrum_str.split()
            spectrum_values = np.array([float(val) for val in spectrum_values])
            return spectrum_values

        # Ask the user to choose the saved session file
        file_path = filedialog.askopenfilename(
            filetypes=[
                ("Shape Files", "*.shp"),  # Filter for Shapefiles
                ("GeoPackage Files", "*.gpkg"),  # Filter for GeoPackage files
            ]
        )

        if file_path:
            if file_path.endswith(".shp"):
                gpf = gpd.read_file(file_path)
            elif file_path.endswith(".gpkg"):
                gpf = gpd.read_file(file_path, driver="GPKG")
                mean_spectrum = gpf["Spectrum Mean"].apply(parse_spectrum).tolist()
                parameter_values = (
                    gpf["Parameters Mean"].apply(parse_parameter_values).tolist()
                )
                # if 'Ratio Spectrum' in gpf.columns:
                #     ratio_spectra = gpf['Ratio Spectrum'].apply(parse_spectrum).tolist()
            else:
                # Handle unsupported file types
                print("Unsupported file type. Please select a .shp or .gpkg file.")
                return
            self.restored_min_id_lines = []
            for i, geom in enumerate(gpf.geometry):
                print(f"reading in polygon {i}")
                # (lon, lat) <--> (x,y)

                # Handle both Polygon and MultiPolygon geometries
                if geom.geom_type == 'MultiPolygon':
                    # For MultiPolygon, process each component polygon
                    for poly in geom.geoms:
                        polygon_coords = poly.exterior.coords.xy
                        polygon_pixel_coords = []
                        # translate from geographic to pixel coordinates
                        for lon, lat in zip(polygon_coords[0], polygon_coords[1]):
                            y, x = self.left_rio.index(lon, lat)
                            polygon_pixel_coords.append((x, y))
                        self.all_polygons.append(polygon_pixel_coords)
                        self.polygon_colors.append(gpf.Color[i])
                        self.restored_min_id_lines.append(
                            [
                                gpf.get("Mineral ID 1")[i],
                                gpf.get("Mineral ID 2")[i],
                                gpf.get("Mineral ID 3")[i],
                                gpf.get("Mineral ID 4")[i],
                            ]
                        )
                        # if the gpf has 'Best Denom ID' column, populate the denominator index
                        if "Best Denom ID" in gpf.columns:
                            self.best_denom_id_list.append(gpf.get("Best Denom ID")[i])
                            print(f'Best Denom ID: {gpf.get("Best Denom ID")[i]}')
                        else:
                            self.best_denom_id_list.append(np.nan)
                        if "Ratio Spectrum" in gpf.columns:
                            self.ratio_spectrum_list.append(gpf.get("Ratio Spectrum")[i])
                        else:
                            self.ratio_spectrum_list.append(np.nan)

                        self.polygon_spectra.append(mean_spectrum[i])
                        self.polygon_params.append(parameter_values[i])
                else:
                    # Handle simple Polygon
                    polygon_coords = geom.exterior.coords.xy
                    polygon_pixel_coords = []
                    # translate from geographic to pixel coordinates
                    for lon, lat in zip(polygon_coords[0], polygon_coords[1]):
                        y, x = self.left_rio.index(lon, lat)
                        polygon_pixel_coords.append((x, y))
                    self.all_polygons.append(polygon_pixel_coords)
                    self.polygon_colors.append(gpf.Color[i])
                    self.restored_min_id_lines.append(
                        [
                            gpf.get("Mineral ID 1")[i],
                            gpf.get("Mineral ID 2")[i],
                            gpf.get("Mineral ID 3")[i],
                            gpf.get("Mineral ID 4")[i],
                        ]
                    )
                    # if the gpf has 'Best Denom ID' column, populate the denominator index
                    if "Best Denom ID" in gpf.columns:
                        self.best_denom_id_list.append(gpf.get("Best Denom ID")[i])
                        print(f'Best Denom ID: {gpf.get("Best Denom ID")[i]}')
                    else:
                        self.best_denom_id_list.append(np.nan)
                    if "Ratio Spectrum" in gpf.columns:
                        self.ratio_spectrum_list.append(gpf.get("Ratio Spectrum")[i])
                    else:
                        self.ratio_spectrum_list.append(np.nan)

                    self.polygon_spectra.append(mean_spectrum[i])
                    self.polygon_params.append(parameter_values[i])

            # self.extract_spectra_from_polygons()
            self.draw_all_polygons()
            print(self.restored_min_id_lines)
            for idx, line in enumerate(self.restored_min_id_lines):
                print(idx)
                # Convert tuple to list, modify, and convert back to tuple
                self.polygons_table_data[idx] = list(self.polygons_table_data[idx])
                self.polygons_table_data[idx][5:9] = line  # Assign values
                self.polygons_table_data[idx] = tuple(self.polygons_table_data[idx])
            self.update_polygons_table()

    # ----------------------------------------------------------------
    # Collecting point spectra
    # ----------------------------------------------------------------
    # todo:
    # 1. fix the y-axis scales when library spectra are plotted
    def update_collected_points_table(self, row=None):

        if hasattr(self, "collected_points_tree"):
            if self.collected_points_window.winfo_exists():
                self.collected_points_tree.insert("", "end", values=row)

            else:
                self.display_collected_points()
        else:
            self.display_collected_points()

    def display_collected_points(self):
        # Create a new window
        self.collected_points_window = tk.Toplevel(self.root)
        self.collected_points_window.title("Collected Points")

        # Create a Treeview widget
        self.collected_points_tree = ttk.Treeview(self.collected_points_window)

        # Get the column names from the DataFrame
        columns = self.collected_points.columns.tolist()

        # Configure the Treeview widget
        self.collected_points_tree["columns"] = columns
        for col in columns:
            self.collected_points_tree.heading(col, text=col)
            self.collected_points_tree.column(col, width=100)

        # Add the data from the DataFrame to the Treeview
        for _, row in self.collected_points.iterrows():
            self.collected_points_tree.insert("", "end", values=row.tolist())

        # Pack the Treeview widget into the window
        self.collected_points_tree.pack()

        # Add a button to save the DataFrame to an Excel file
        save_button = tk.Button(
            self.collected_points_window,
            text="Save to Excel",
            command=self.save_collected_points_to_excel,
        )
        save_button.pack()

        # Add a button to plot the collected spectra
        plot_button = tk.Button(
            self.collected_points_window,
            text="Plot Spectra",
            command=self.plot_collected_spectra,
        )
        plot_button.pack()

        self.collected_points_tree.bind("<Double-1>", self.edit_points_tree_item)
        self.collected_points_tree.bind(
            "<Button-2>", self.show_collected_points_tree_context_menu
        )

        self.collected_points_tree_context_menu = tk.Menu(
            self.collected_points_tree, tearoff=0
        )
        self.collected_points_tree_context_menu.add_command(
            label="Delete point", command=self.delete_collected_point
        )

        # Show the window
        # self.collected_points_window.mainloop()

    def edit_points_tree_item(self, event):
        # Get the selected item and column
        selected_item = self.collected_points_tree.identify_row(event.y)
        column = self.collected_points_tree.identify_column(event.x)

        # Check if the column is editable
        if column in ["#2", "#3", "#4", "#5", "#6"]:
            # Create an Entry widget
            self.entry = tk.Entry(self.collected_points_tree)
            self.entry.insert(
                0,
                self.collected_points_tree.item(selected_item, "values")[
                    int(column[1:]) - 1
                ],
            )

            # Place the Entry widget over the cell
            self.entry.place(
                x=event.x,
                y=event.y,
                width=self.collected_points_tree.column(column, "width"),
            )

            # Bind the Return key and the focus out event to the save_edit method
            self.entry.bind(
                "<Return>",
                lambda _: self.save_collected_point_table_edit(selected_item, column),
            )
            self.entry.bind(
                "<FocusOut>",
                lambda _: self.save_collected_point_table_edit(selected_item, column),
            )

    def save_collected_point_table_edit(self, item, column):
        # Get the new value
        new_value = self.entry.get()

        # Update the Treeview
        self.collected_points_tree.set(item, column, new_value)

        # Update the DataFrame
        self.collected_points.loc[
            self.collected_points["name"]
            == self.collected_points_tree.item(item, "values")[0],
            self.collected_points_tree.column(column, "id"),
        ] = new_value

        # Destroy the Entry widget
        self.entry.destroy()

    def show_collected_points_tree_context_menu(self, event):
        self.tmp_event = event
        item = self.collected_points_tree.identify_row(event.y)
        if item:
            self.collected_points_tree.selection_set(item)
            self.collected_points_tree_context_menu.post(event.x_root, event.y_root)

    def save_collected_points_to_excel(self):
        # Ask the user to choose the file name and location
        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")]
        )

        # Save the DataFrame to an Excel file
        if file_path:
            self.collected_points.to_excel(file_path)

    def delete_collected_point(self):
        selected_item = self.collected_points_tree.selection()[0]  # get selected item
        item_values = self.collected_points_tree.item(selected_item, "values")
        point_name = item_values[0]  # assuming the first column is the point name
        self.collected_points = self.collected_points[
            self.collected_points.name != point_name
        ]
        self.collected_points_tree.delete(selected_item)

    def collect_point(self):
        column_headers = self.collected_points.columns.tolist()
        point_index = len(self.collected_points)
        point_name = f"{self.img_base_name}_{point_index}"
        x, y = self.current_spec_loc
        spec = (self.left_wvl, self.spectrum)
        params = (self.right_wvl, self.parameter_values)
        # get geolocation from x, y
        lon, lat = self.left_rio.xy(y, x)
        if self.denom_spectrum is not None:
            denom = self.denom_spectrum
        else:
            denom = ""
        row = (point_name, "", "", "", "", "", x, y, lon, lat, spec, denom, params)
        row_dict = {}
        for key, value in zip(column_headers, row):
            row_dict[key] = value
        # add row to the self.collected_points data frame
        self.collected_points.loc[len(self.collected_points)] = row_dict

        self.update_collected_points_table(row=row)

    def plot_collected_spectra(self):
        # Create a new window
        self.collected_spectral_plot_window = tk.Toplevel(self.root)
        self.collected_spectral_plot_window.title("Collected Spectra")

        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.collected_spectral_plot_window)
        ui_frame.pack(side=tk.RIGHT, fill=tk.BOTH)

        # Create a Figure and a FigureCanvasTkAgg object
        fig = Figure(figsize=(5, 4), dpi=100)
        self.collected_spectra_canvas = FigureCanvasTkAgg(
            fig, master=self.collected_spectral_plot_window
        )

        # Create an Axes object and plot the spectra
        ax = fig.add_subplot(111)
        for _, row in self.collected_points.iterrows():
            # retrieve the values from the "spectrum" column
            spec_label = row["name"]
            wvl, spectrum = row["spectrum"]
            color = row["color"]
            if color == "":
                color = "k"
            ax.plot(wvl, spectrum, label=spec_label, color=color)

        # add a button to turn the legend on or off
        self.toggle_collected_legend_button = tk.Button(
            ui_frame,
            text="Legend (On)",
            command=self.toggle_collected_spectral_plot_legend,
        )
        self.toggle_collected_legend_button.pack(side=tk.TOP)

        # Add a button to reset x-axis span
        self.reset_spectral_x_axis_button = tk.Button(
            ui_frame, text="Reset X-Axis Span", command=self.reset_collected_x_axis_span
        )
        self.reset_spectral_x_axis_button.pack(side=tk.TOP)

        # Add built-in span options to a dropdown menu
        span_options = [
            "Full Span",
            "410 - 1000 nm",
            "410 - 2500 nm",
            "1000 - 2600 nm",
            "1200 - 2000 nm",
            "1800 - 2500 nm",
            "2000 - 2500 nm",
            "2700 - 3900 nm",
        ]
        self.collected_span_var = tk.StringVar()
        self.collected_span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(
            ui_frame,
            textvariable=self.collected_span_var,
            values=span_options,
            state="readonly",
        )
        span_menu.pack(side=tk.TOP)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_collected_x_axis_span)

        # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        self.collected_library_spectra_var = tk.StringVar(ui_frame)
        self.collected_library_spectra_var.set("None")  # default value
        self.collected_library_spectra_list = [
            name.split("/")[-1] for name in self.usgs_spectra_folders
        ]
        # sort library_spectra_list alphabetically
        self.collected_library_spectra_list.sort()
        # label the library spectra drop down menu "library spectra"
        self.collected_library_spectra_label = tk.Label(
            ui_frame, text="Library Spectra"
        )
        self.collected_library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.collected_library_spectra_dropdown = tk.OptionMenu(
            ui_frame,
            self.collected_library_spectra_var,
            *self.collected_library_spectra_list,
            command=self.plot_collected_library_spectra,
        )
        self.collected_library_spectra_dropdown.pack(side=tk.TOP)

        # add a button to remove library spectra from the plot
        self.remove_collected_library_spectra_button = tk.Button(
            ui_frame,
            text="Remove Library Spectra",
            command=self.remove_collected_library_spectra,
        )
        self.remove_collected_library_spectra_button.pack(side=tk.TOP)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(
            self.collected_spectra_canvas, self.collected_spectral_plot_window
        )
        toolbar.update()
        self.collected_spectra_canvas.get_tk_widget().pack(
            side=tk.TOP, fill=tk.BOTH, expand=True
        )

        # Add span selector for x-axis
        self.create_collected_x_axis_span_selector()

        # Pack the FigureCanvasTkAgg object into the window
        self.collected_spectra_canvas.draw()
        self.collected_spectra_canvas.get_tk_widget().pack(
            side=tk.TOP, fill=tk.BOTH, expand=1
        )

    def toggle_collected_spectral_plot_legend(self):
        ax = self.collected_spectra_canvas.figure.axes[0]
        if ax.get_legend() is None:
            ax.legend()
            self.toggle_collected_legend_button.config(
                text="Legend (On)", relief="raised"
            )
        else:
            ax.legend_.remove()
            self.toggle_collected_legend_button.config(
                text="Legend (Off)", relief="sunken"
            )
        self.collected_spectra_canvas.draw()

    def reset_collected_x_axis_span(self):
        # Reset x-axis span to the default range
        ax = self.collected_spectra_canvas.figure.axes[0]
        ax.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        ymin = []
        ymax = []
        for line in ax.lines:
            ymin.append(np.nanmin(line.get_ydata()[xmin_idx:xmax_idx]))
            ymax.append(np.nanmax(line.get_ydata()[xmin_idx:xmax_idx]))
        y_min = np.nanmin(ymin)
        y_max = np.nanmax(ymax)
        buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
        ax.set_ylim(y_min - buffer, y_max + buffer)
        self.collected_spectra_canvas.draw()

    def create_collected_x_axis_span_selector(self):
        ax = self.collected_spectra_canvas.figure.axes[0]

        def on_x_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(self.left_wvl - xmin))
            xmax_idx = np.argmin(np.abs(self.left_wvl - xmax))
            ax.set_xlim(self.left_wvl[xmin_idx], self.left_wvl[xmax_idx])

            # Calculate y-limits based on the data within the new span
            ymin = []
            ymax = []
            for line in ax.lines:
                ymin.append(np.nanmin(line.get_ydata()[xmin_idx:xmax_idx]))
                ymax.append(np.nanmax(line.get_ydata()[xmin_idx:xmax_idx]))
            y_min = np.nanmin(ymin)
            y_max = np.nanmax(ymax)

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            ax.set_ylim(y_min - buffer, y_max + buffer)

            self.collected_spectra_canvas.draw()

        self.collected_x_span_selector = SpanSelector(
            ax, on_x_span_select, "horizontal", useblit=True
        )

    def update_collected_x_axis_span(self, event):
        selected_span = self.collected_span_var.get()
        ax = self.collected_spectra_canvas.figure.axes[0]
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_wvl[0], self.left_wvl[-1]),
            "410 - 1000 nm": (410, 1000),
            "410 - 2500 nm": (410, 2500),
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900),
        }

        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            ax.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))

            # Calculate y-limits for each line in the plot based on the data within the new span
            ymin = []
            ymax = []
            for line in ax.lines:
                ymin.append(np.nanmin(line.get_ydata()[xmin_idx:xmax_idx]))
                ymax.append(np.nanmax(line.get_ydata()[xmin_idx:xmax_idx]))

            y_min = np.nanmin(ymin)
            y_max = np.nanmax(ymax)

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            ax.set_ylim(y_min - buffer, y_max + buffer)

            self.collected_spectra_canvas.draw()

    def plot_collected_library_spectra(self, event):
        if self.collected_library_spectra_var.get() != "None":
            # plot the selected library spectrum
            ax = self.collected_spectra_canvas.figure.axes[0]
            path_to_spec_data = self.usgs_spectra_path + event + "/" + event + ".txt"
            path_to_wvl_data = (
                self.usgs_spectra_path + event + "/" + "*Wavelengths*" + ".txt"
            )
            library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
            library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            # plot the values in spectrum_df on self.spectral_ax
            ax.plot(
                library_wvl,
                library_reflectance,
                label=self.collected_library_spectra_var.get(),
            )
            # scale each line so that the min and max values equal the min and max max value of self.spectrum
            for line in ax.lines[1:]:
                new_y_data = line.get_ydata() * (
                    np.nanmax(self.spectrum) / np.nanmax(line.get_ydata())
                    - np.nanmin(self.spectrum)
                ) + np.nanmin(self.spectrum)
                line.set_ydata(new_y_data)
            # min_y, max_y = np.nanmin(library_reflectance), np.nanmax(library_reflectance)
            # buffer = (max_y - min_y) * 0.1
            # new_y_lim = (np.nanmin((ymin, min_y - buffer)), np.nanmax((ymax, max_y + buffer)))
            ax.set_ylim((ymin, ymax))
            ax.set_xlim(xmin, xmax)
            ax.legend()
            self.collected_spectra_canvas.draw()

    def remove_collected_library_spectra(self):
        # remove the selected library spectrum
        ax = self.collected_spectra_canvas.figure.axes[0]
        ax.lines[-1].remove()
        if ax.legend_ is not None:
            ax.legend_.remove()
        ax.legend()
        self.collected_spectra_canvas.draw()

    # ----------------------------------------------------------------
    # Spectral plotting area
    # ----------------------------------------------------------------
    # polygons
    # ----------------------------------------------------------------
    def create_polygons_spectral_plot(self):

        if not self.polygon_spectra:
            self.extract_spectra_from_polygons()
        else:
            self.polygons_spectral_window = tk.Toplevel(self.root)
            self.polygons_spectral_window.title("ROI Spectral Plot")

            # Create a frame to hold UI elements with a fixed size
            ui_frame = tk.Frame(self.polygons_spectral_window)
            ui_frame.pack(side=tk.RIGHT, fill=tk.X)

            polygons_spectral_figure, self.polygons_spectral_ax = plt.subplots(
                figsize=(5, 3)
            )

            for poly_color, s in zip(self.polygon_colors, self.polygon_spectra):
                s = np.array(s)
                s = s.flatten()
                self.polygon_spectral_lines.append(
                    self.polygons_spectral_ax.plot(self.left_wvl, s, color=poly_color)
                )
            # self.polygon_spectral_lines, = self.polygons_spectral_ax.plot(self.left_wvl, self.polygon_spectra)
            xmin, xmax = self.polygons_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            all_min_y = []
            all_max_y = []
            for s in self.polygon_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                    s[xmin_idx:xmax_idx]
                )
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_spectral_canvas = FigureCanvasTkAgg(
                polygons_spectral_figure, master=self.polygons_spectral_window
            )
            self.polygons_spectral_canvas.get_tk_widget().pack(
                fill=tk.BOTH, expand=True
            )

            # Add a button to reset x-axis span
            self.reset_polygons_x_axis_button = tk.Button(
                ui_frame,
                text="Reset X-Axis Span",
                command=self.reset_polygons_x_axis_span,
            )
            self.reset_polygons_x_axis_button.pack(side=tk.TOP)

            # Add a button to ratio spectra to a a chosen denominator spectrum
            self.ratio_spectra_button = tk.Button(
                ui_frame, text="Ratio Spectra", command=self.ratio_polygon_spectra
            )
            self.ratio_spectra_button.pack(side=tk.TOP)

            # Add built-in span options to a dropdown menu
            span_options = [
                "Full Span",
                "410 - 1000 nm",
                "410 - 2500 nm",
                "1000 - 2600 nm",
                "1200 - 2000 nm",
                "1800 - 2500 nm",
                "2000 - 2500 nm",
                "2700 - 3900 nm",
            ]
            self.polygons_span_var = tk.StringVar()
            self.polygons_span_var.set("Full Span")  # Set the default span option
            span_menu = ttk.Combobox(
                ui_frame,
                textvariable=self.polygons_span_var,
                values=span_options,
                state="readonly",
            )
            span_menu.pack(side=tk.RIGHT)

            # add a button to turn the legend on or off
            self.toggle_polygons_legend_button = tk.Button(
                ui_frame,
                text="Legend (On)",
                command=self.toggle_polygons_spectral_plot_legend,
            )
            self.toggle_polygons_legend_button.pack(side=tk.TOP)

            # add a button to ignore bad bands
            self.polygons_ignore_bad_bands_button = tk.Button(
                ui_frame,
                text="Ignore Bad Bands (Off)",
                command=self.toggle_polygons_ignore_bad_bands,
            )
            self.polygons_ignore_bad_bands_button.pack(side=tk.TOP)

            # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
            self.polygons_library_spectra_var = tk.StringVar(ui_frame)
            self.polygons_library_spectra_var.set("None")  # default value
            self.polygons_library_spectra_list = [
                name.split("/")[-1] for name in self.usgs_spectra_folders
            ]
            # sort library_spectra_list alphabetically
            self.polygons_library_spectra_list.sort()
            # label the library spectra drop down menu "library spectra"
            self.polygons_library_spectra_label = tk.Label(
                ui_frame, text="Library Spectra"
            )
            self.polygons_library_spectra_label.pack(side=tk.TOP)
            # create the drop down menu
            self.polygons_library_spectra_dropdown = tk.OptionMenu(
                ui_frame,
                self.polygons_library_spectra_var,
                *self.polygons_library_spectra_list,
                command=self.polygons_plot_library_spectra,
            )
            self.polygons_library_spectra_dropdown.pack(side=tk.TOP)

            # add a button to remove library spectra from the plot
            self.polygons_remove_library_spectra_button = tk.Button(
                ui_frame,
                text="Remove Library Spectra",
                command=self.polygons_remove_library_spectra,
            )
            self.polygons_remove_library_spectra_button.pack(side=tk.TOP)

            # Bind an event to update the x-axis span when a span option is selected
            span_menu.bind("<<ComboboxSelected>>", self.update_polygons_x_axis_span)

            # Create a toolbar for the spectral plot
            toolbar = NavigationToolbar2Tk(
                self.polygons_spectral_canvas, self.polygons_spectral_window
            )
            toolbar.update()
            self.polygons_spectral_canvas.get_tk_widget().pack(
                side=tk.TOP, fill=tk.BOTH, expand=True
            )

            # Add span selector for x-axis
            self.create_polygons_x_axis_span_selector(self.left_wvl)

    def update_polygons_spectral_plot(self, use_unified: bool = True):
        """
        Update the ROI spectral plot.

        Args:
            use_unified: If True, use the unified plotting system (SpectralPlotWindow).
                        If False, use the legacy separate window.
        """
        if use_unified and self.polygon_spectra:
            self._add_roi_spectra_to_active_plot()
            return

        # Legacy behavior below
        if (
            self.polygons_spectral_window is None
            or not self.polygons_spectral_window.winfo_exists()
        ):
            self.polygon_spectral_lines = []
            self.create_polygons_spectral_plot()
        elif self.polygon_spectra:
            # clear everything from the plot before plotting
            self.polygons_spectral_ax.clear()
            self.polygon_spectral_lines = []

            for c, s in zip(self.polygon_colors, self.polygon_spectra):
                s = np.array(s)
                s = s.flatten()
                if self.ignore_polygons_bad_bands_flag:
                    s = np.where(np.array(self.bad_bands) == 0, np.nan, s)
                self.polygon_spectral_lines.append(
                    self.polygons_spectral_ax.plot(self.left_wvl, s, color=c)
                )
            xmin, xmax = self.polygons_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            all_min_y = []
            all_max_y = []
            for s in self.polygon_spectra:
                # Check if the spectrum slice has any valid data
                spectrum_slice = s[xmin_idx:xmax_idx]
                if not np.all(np.isnan(spectrum_slice)):
                    min_y, max_y = np.nanmin(spectrum_slice), np.nanmax(spectrum_slice)
                    all_min_y.append(min_y)
                    all_max_y.append(max_y)

            # Only set limits if we have valid data
            if all_min_y and all_max_y:
                min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
                buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
                self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.polygons_spectral_canvas.draw()
        else:
            self.polygons_spectral_ax.clear()
            self.polygons_spectral_canvas.draw()

    def toggle_polygons_spectral_plot_legend(self):
        ax = self.polygons_spectral_canvas.figure.axes[0]
        if ax.get_legend() is None:
            ax.legend()
            self.toggle_polygons_legend_button.config(
                text="Legend (On)", relief="raised"
            )
        else:
            ax.legend_.remove()
            self.toggle_polygons_legend_button.config(
                text="Legend (Off)", relief="sunken"
            )
        self.polygons_spectral_canvas.draw()

    def toggle_polygons_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_polygons_bad_bands_flag:
                self.ignore_polygons_bad_bands_flag = False
                self.polygons_ignore_bad_bands_button.config(
                    text="Ignore Bad Bands (Off)", relief="sunken"
                )
            else:
                self.ignore_polygons_bad_bands_flag = True
                self.polygons_ignore_bad_bands_button.config(
                    text="Ignore Bad Bands (On)", relief="raised"
                )
            self.update_polygons_spectral_plot()

    def polygons_plot_library_spectra(self, event):
        if self.polygons_library_spectra_var.get() != "None":
            # plot the selected library spectrum
            ax = self.polygons_spectral_canvas.figure.axes[0]
            path_to_spec_data = (
                self.usgs_spectra_path
                + self.polygons_library_spectra_var.get()
                + "/"
                + self.polygons_library_spectra_var.get()
                + ".txt"
            )
            path_to_wvl_data = (
                self.usgs_spectra_path
                + self.polygons_library_spectra_var.get()
                + "/"
                + "*Wavelengths*"
                + ".txt"
            )
            library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
            library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            # plot the values in spectrum_df on self.spectral_ax
            ax.plot(
                library_wvl,
                library_reflectance,
                label=self.polygons_library_spectra_var.get(),
            )
            self.polygons_library_reflectance.append(ax.lines[-1])
            # scale each line so that the min and max values equal the min and max max value of self.spectrum
            ymin_ = []
            ymax_ = []
            for line in self.polygon_spectral_lines:
                ymin_.append(np.nanmin(line[0].get_ydata()[xmin_idx:xmax_idx]))
                ymax_.append(np.nanmax(line[0].get_ydata()[xmin_idx:xmax_idx]))
            y_min = np.nanmin(ymin_)
            y_max = np.nanmax(ymax_)

            for line in self.polygons_library_reflectance:
                new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * (
                    (y_max - y_min)
                    / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))
                ) + y_min
                line.set_ydata(new_y_data)
            # min_y, max_y = np.nanmin(library_reflectance), np.nanmax(library_reflectance)
            # buffer = (max_y - min_y) * 0.1
            # new_y_lim = (np.nanmin((ymin, min_y - buffer)), np.nanmax((ymax, max_y + buffer)))
            ax.set_ylim((ymin, ymax))
            ax.set_xlim(xmin, xmax)
            ax.legend()
            self.polygons_spectral_canvas.draw()

    def polygons_remove_library_spectra(self):
        # remove the selected library spectrum
        ax = self.polygons_spectral_canvas.figure.axes[0]
        self.polygons_library_reflectance[-1].remove()
        if ax.legend_ is not None:
            ax.legend_.remove()
        ax.legend()
        self.polygons_library_reflectance.pop()
        self.polygons_spectral_canvas.draw()

    def reset_polygons_x_axis_span(self):
        # Reset x-axis span to the default range
        self.polygons_spectral_ax.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.polygons_spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        all_min_y = []
        all_max_y = []
        for s in self.polygon_spectra:
            min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                s[xmin_idx:xmax_idx]
            )
            all_min_y.append(min_y)
            all_max_y.append(max_y)
        min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
        self.polygons_spectral_canvas.draw()

    def update_polygons_x_axis_span(self, event):
        selected_span = self.polygons_span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_wvl[0], self.left_wvl[-1]),
            "410 - 1000 nm": (410, 1000),
            "410 - 2500 nm": (410, 2500),
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900),
        }
        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            self.polygons_spectral_ax.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))

            # Calculate y-limits based on the data within the new span
            all_min_y = []
            all_max_y = []
            for s in self.polygon_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                    s[xmin_idx:xmax_idx]
                )
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)

            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_spectral_canvas.draw()

    def create_polygons_x_axis_span_selector(self, x_data):
        def on_x_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(x_data - xmin))
            xmax_idx = np.argmin(np.abs(x_data - xmax))
            self.polygons_spectral_ax.set_xlim(x_data[xmin_idx], x_data[xmax_idx])

            # Calculate y-limits based on the data within the new span
            all_min_y = []
            all_max_y = []
            for s in self.polygon_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                    s[xmin_idx:xmax_idx]
                )
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)

            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_spectral_canvas.draw()

        self.x_polygons_span_selector = SpanSelector(
            self.polygons_spectral_ax, on_x_span_select, "horizontal", useblit=True
        )

        # ----------------------------------------------------------------
        # Polygon ratios
        # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    #  Polygon ratio spectra
    # ----------------------------------------------------------------
    def ratio_polygon_spectra(self):
        if self.denominator_index < 0 or self.denominator_index >= len(
            self.polygon_spectra
        ):
            print("Invalid denominator index")
            return

        denominator_spectrum = self.polygon_spectra[self.denominator_index]
        self.polygon_ratio_spectra = []
        for spectrum in self.polygon_spectra:
            self.polygon_ratio_spectra.append(spectrum / denominator_spectrum)

        # Update the spectral plot
        self.update_polygons_ratio_spectral_plot()

    def create_polygons_ratio_spectral_plot(self):

        if not self.polygon_ratio_spectra:
            self.extract_spectra_from_polygons()
            self.ratio_polygon_spectra()
        else:
            self.polygons_ratio_spectral_window = tk.Toplevel(self.root)
            self.polygons_ratio_spectral_window.title("ROI Spectral Ratio Plot")

            # Create a frame to hold UI elements with a fixed size
            ui_frame = tk.Frame(self.polygons_ratio_spectral_window)
            ui_frame.pack(side=tk.RIGHT, fill=tk.X)

            polygons_ratio_spectral_figure, self.polygons_ratio_spectral_ax = (
                plt.subplots(figsize=(5, 3))
            )

            for i, (poly_color, s) in enumerate(
                zip(self.polygon_colors, self.polygon_ratio_spectra)
            ):
                if i == self.denominator_index:
                    continue
                s = s.flatten()
                l = self.polygons_ratio_spectral_ax.plot(
                    self.left_wvl, s, color=poly_color
                )
                self.polygon_ratio_spectral_lines.append(l)
            # self.polygon_ratio_spectral_lines, = self.polygons_spectral_ax.plot(self.left_wvl, self.polygon_spectra)
            xmin, xmax = self.polygons_ratio_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            all_min_y = []
            all_max_y = []
            for s in self.polygon_ratio_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                    s[xmin_idx:xmax_idx]
                )
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_ratio_spectral_canvas = FigureCanvasTkAgg(
                polygons_ratio_spectral_figure,
                master=self.polygons_ratio_spectral_window,
            )
            self.polygons_ratio_spectral_canvas.get_tk_widget().pack(
                fill=tk.BOTH, expand=True
            )

            # add a button to turn the legend on or off
            self.toggle_polygons_ratio_legend_button = tk.Button(
                ui_frame,
                text="Legend (On)",
                command=self.toggle_polygons_ratio_spectral_plot_legend,
            )
            self.toggle_polygons_ratio_legend_button.pack(side=tk.TOP)

            # Add a button to reset x-axis span
            self.reset_polygons_ratio_x_axis_button = tk.Button(
                ui_frame,
                text="Reset X-Axis Span",
                command=self.reset_polygons_ratio_x_axis_span,
            )
            self.reset_polygons_ratio_x_axis_button.pack(side=tk.TOP)

            # Add built-in span options to a dropdown menu
            span_options = [
                "Full Span",
                "410 - 1000 nm",
                "410 - 2500 nm",
                "1000 - 2600 nm",
                "1200 - 2000 nm",
                "1800 - 2500 nm",
                "2000 - 2500 nm",
                "2700 - 3900 nm",
            ]
            self.polygons_ratio_span_var = tk.StringVar()
            self.polygons_ratio_span_var.set("Full Span")  # Set the default span option
            span_menu = ttk.Combobox(
                ui_frame,
                textvariable=self.polygons_ratio_span_var,
                values=span_options,
                state="readonly",
            )
            span_menu.pack(side=tk.TOP)

            # Bind an event to update the x-axis span when a span option is selected
            span_menu.bind(
                "<<ComboboxSelected>>", self.update_polygons_ratio_x_axis_span
            )

            # add a button to ignore bad bands
            self.polygons_ratio_ignore_bad_bands_button = tk.Button(
                ui_frame,
                text="Ignore Bad Bands (Off)",
                command=self.toggle_polygons_ratio_ignore_bad_bands,
            )
            self.polygons_ratio_ignore_bad_bands_button.pack(side=tk.TOP)

            # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
            self.polygons_ratio_library_spectra_var = tk.StringVar(ui_frame)
            self.polygons_ratio_library_spectra_var.set("None")  # default value
            self.polygons_ratio_library_spectra_list = [
                name.split("/")[-1] for name in self.usgs_spectra_folders
            ]
            # sort library_spectra_list alphabetically
            self.polygons_ratio_library_spectra_list.sort()
            # label the library spectra drop down menu "library spectra"
            self.polygons_ratio_library_spectra_label = tk.Label(
                ui_frame, text="Library Spectra"
            )
            self.polygons_ratio_library_spectra_label.pack(side=tk.TOP)
            # create the drop down menu
            self.polygons_ratio_library_spectra_dropdown = tk.OptionMenu(
                ui_frame,
                self.polygons_ratio_library_spectra_var,
                *self.polygons_ratio_library_spectra_list,
                command=self.polygons_ratio_plot_library_spectra,
            )
            self.polygons_ratio_library_spectra_dropdown.pack(side=tk.TOP)

            # add a button to remove library spectra from the plot
            self.polygons_ratio_remove_library_spectra_button = tk.Button(
                ui_frame,
                text="Remove Library Spectra",
                command=self.polygons_ratio_remove_library_spectra,
            )
            self.polygons_ratio_remove_library_spectra_button.pack(side=tk.TOP)

            # Create a toolbar for the spectral plot
            toolbar = NavigationToolbar2Tk(
                self.polygons_ratio_spectral_canvas, self.polygons_ratio_spectral_window
            )
            toolbar.update()
            self.polygons_ratio_spectral_canvas.get_tk_widget().pack(
                side=tk.TOP, fill=tk.BOTH, expand=True
            )

            # Add span selector for x-axis
            self.create_polygons_ratio_x_axis_span_selector(self.left_wvl)

    def toggle_polygons_ratio_spectral_plot_legend(self):
        ax = self.polygons_ratio_spectral_canvas.figure.axes[0]
        if ax.get_legend() is None:
            ax.legend()
            self.toggle_polygons_ratio_legend_button.config(
                text="Legend (On)", relief="raised"
            )
        else:
            ax.legend_.remove()
            self.toggle_polygons_ratio_legend_button.config(
                text="Legend (Off)", relief="sunken"
            )
        self.polygons_ratio_spectral_canvas.draw()

    def toggle_polygons_ratio_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_bad_bands_flag:
                self.ignore_bad_bands_flag = False
                self.polygons_ratio_ignore_bad_bands_button.config(
                    text="Ignore Bad Bands (Off)", relief="sunken"
                )
            else:
                self.ignore_bad_bands_flag = True
                self.polygons_ratio_ignore_bad_bands_button.config(
                    text="Ignore Bad Bands (On)", relief="raised"
                )
            self.update_polygons_ratio_spectral_plot()

    def polygons_ratio_plot_library_spectra(self, event):
        if self.polygons_ratio_library_spectra_var.get() != "None":
            # plot the selected library spectrum
            ax = self.polygons_ratio_spectral_canvas.figure.axes[0]
            path_to_spec_data = (
                self.usgs_spectra_path
                + self.polygons_ratio_library_spectra_var.get()
                + "/"
                + self.polygons_ratio_library_spectra_var.get()
                + ".txt"
            )
            path_to_wvl_data = (
                self.usgs_spectra_path
                + self.polygons_ratio_library_spectra_var.get()
                + "/"
                + "*Wavelengths*"
                + ".txt"
            )
            library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
            library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            # plot the values in spectrum_df on self.spectral_ax
            ax.plot(
                library_wvl,
                library_reflectance,
                label=self.polygons_ratio_library_spectra_var.get(),
            )
            self.polygons_ratio_library_reflectance.append(ax.lines[-1])
            # scale each line so that the min and max values equal the min and max max value of self.spectrum
            ymin_ = []
            ymax_ = []
            for line in self.polygon_ratio_spectral_lines:
                ymin_.append(np.nanmin(line[0].get_ydata()[xmin_idx:xmax_idx]))
                ymax_.append(np.nanmax(line[0].get_ydata()[xmin_idx:xmax_idx]))
            y_min = np.nanmin(ymin_)
            y_max = np.nanmax(ymax_)

            for line in self.polygons_ratio_library_reflectance:
                new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * (
                    (y_max - y_min)
                    / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))
                ) + y_min
                line.set_ydata(new_y_data)
            # min_y, max_y = np.nanmin(library_reflectance), np.nanmax(library_reflectance)
            # buffer = (max_y - min_y) * 0.1
            # new_y_lim = (np.nanmin((ymin, min_y - buffer)), np.nanmax((ymax, max_y + buffer)))
            ax.set_ylim((ymin, ymax))
            ax.set_xlim(xmin, xmax)
            ax.legend()
            self.polygons_ratio_spectral_canvas.draw()

    def polygons_ratio_remove_library_spectra(self):
        # remove the selected library spectrum
        ax = self.polygons_ratio_spectral_canvas.figure.axes[0]
        self.polygons_ratio_library_reflectance[-1].remove()
        if ax.legend_ is not None:
            ax.legend_.remove()
        ax.legend()
        self.polygons_ratio_library_reflectance.pop()
        self.polygons_ratio_spectral_canvas.draw()

    def update_polygons_ratio_spectral_plot(self):
        if (
            self.polygons_ratio_spectral_window is None
            or not self.polygons_ratio_spectral_window.winfo_exists()
        ):
            self.polygon_ratio_spectral_lines = []
            self.create_polygons_ratio_spectral_plot()
        elif self.polygon_ratio_spectra:
            # clear everything from the plot before plotting
            self.polygons_ratio_spectral_ax.clear()
            self.polygon_ratio_spectral_lines = []

            for i, (c, s) in enumerate(
                zip(self.polygon_colors, self.polygon_ratio_spectra)
            ):
                if i == self.denominator_index:
                    continue
                s = s.flatten()
                if self.ignore_bad_bands_flag:
                    s = np.where(np.array(self.bad_bands) == 0, np.nan, s)
                l = self.polygons_ratio_spectral_ax.plot(self.left_wvl, s, color=c)
                self.polygon_ratio_spectral_lines.append(l)
            xmin, xmax = self.polygons_ratio_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))

            all_min_y = []
            all_max_y = []
            for s in self.polygon_ratio_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                    s[xmin_idx:xmax_idx]
                )
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.polygons_ratio_spectral_canvas.draw()
        else:
            self.polygons_ratio_spectral_ax.clear()
            self.polygons_ratio_spectral_canvas.draw()

    def reset_polygons_ratio_x_axis_span(self):
        # Reset x-axis span to the default range
        self.polygons_ratio_spectral_ax.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.polygons_ratio_spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        all_min_y = []
        all_max_y = []
        for i, s in enumerate(self.polygon_ratio_spectra):
            if i == self.denominator_index:
                continue
            min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                s[xmin_idx:xmax_idx]
            )
            all_min_y.append(min_y)
            all_max_y.append(max_y)
        min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.polygons_ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
        self.polygons_ratio_spectral_canvas.draw()

    def update_polygons_ratio_x_axis_span(self, event):
        selected_span = self.polygons_ratio_span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_wvl[0], self.left_wvl[-1]),
            "410 - 1000 nm": (410, 1000),
            "410 - 2500 nm": (410, 2500),
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900),
        }
        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            self.polygons_ratio_spectral_ax.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))

            # Calculate y-limits based on the data within the new span
            all_min_y = []
            all_max_y = []
            for s in self.polygon_ratio_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                    s[xmin_idx:xmax_idx]
                )
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)

            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_ratio_spectral_canvas.draw()

    def create_polygons_ratio_x_axis_span_selector(self, x_data):
        def on_x_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(x_data - xmin))
            xmax_idx = np.argmin(np.abs(x_data - xmax))
            self.polygons_ratio_spectral_ax.set_xlim(x_data[xmin_idx], x_data[xmax_idx])

            # Calculate y-limits based on the data within the new span
            all_min_y = []
            all_max_y = []
            for i, s in enumerate(self.polygon_ratio_spectra):
                if i == self.denominator_index:
                    continue
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(
                    s[xmin_idx:xmax_idx]
                )
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)

            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_ratio_spectral_canvas.draw()

        self.x_polygons_ratio_span_selector = SpanSelector(
            self.polygons_ratio_spectral_ax,
            on_x_span_select,
            "horizontal",
            useblit=True,
        )

    # ----------------------------------------------------------------
    # points
    # ----------------------------------------------------------------
    def create_new_spectral_plot_window(self):
        """Create a new empty spectral plot window (menu action)."""
        plot_window = self._create_spectral_plot_window()
        plot_window.focus()
        return plot_window

    def _create_spectral_plot_window(self, title: str = "Spectral Plot") -> SpectralPlotWindow:
        """
        Create and register a new SpectralPlotWindow.

        Returns:
            The created SpectralPlotWindow instance
        """
        def on_close(plot):
            self.plot_manager.unregister_plot(plot.plot_id)

        def on_focus(plot):
            self.plot_manager.set_active_plot(plot.plot_id)

        plot_window = SpectralPlotWindow(
            root=self.root,
            title=title,
            on_close_callback=on_close,
            on_focus_callback=on_focus,
            # Use a provider, not a snapshot: select_bad_bands /
            # reset_bad_bands / on_bad_band_change rebind self.bad_bands
            # to a fresh list, which would leave a snapshot stale and
            # the toggle silently inert.
            bad_bands_provider=lambda: self.bad_bands,
            library_spectra_folders=self.usgs_spectra_folders,
            mica_spectra_names=self.mica_spectra_names,
            usgs_spectra_path=self.usgs_spectra_path,
            mica_spectra_path=self.mica_spectra_path,
        )

        plot_id = self.plot_manager.register_plot(plot_window)
        plot_window.set_plot_id(plot_id)

        return plot_window

    def create_spectral_plot(self):
        """Create a new spectral plot window using the new system."""
        # Create a new SpectralPlotWindow
        plot_window = self._create_spectral_plot_window()

        # Add the current spectrum to it
        spectrum_data = SpectrumData(
            wavelengths=np.array(self.left_wvl),
            values=self.spectrum.copy(),
            label=self.spectrum_label,
        )
        plot_window.add_spectrum(spectrum_data)

        # Store reference for backwards compatibility with old code
        self.spectral_window = plot_window.window
        self.spectral_ax = plot_window._ax
        self.spectral_canvas = plot_window._canvas
        # Store the line for update_spectral_plot compatibility
        if plot_window._spectrum_lines:
            self.spectral_line = plot_window._spectrum_lines[0]

        # Also keep reference to the plot window object
        self._current_point_plot_window = plot_window

        # Legacy: Initialize offset correction flag
        self.offset_correction_applied = False
        self.original_spectrum = self.spectrum.copy()

        # # Create a new window for the spectral plot
        # self.spectral_window = tk.Toplevel(self.root)
        # self.spectral_window.title("Spectral Plot")

        # # Create a frame to hold UI elements with a fixed size
        # ui_frame = tk.Frame(self.spectral_window)
        # ui_frame.pack(side=tk.RIGHT, fill=tk.BOTH)

        # spectral_figure, self.spectral_ax = plt.subplots(figsize=(5,3))
        # self.spectral_line, = self.spectral_ax.plot(self.left_wvl, self.spectrum, label=self.spectrum_label)
        # self.spectral_ax.legend(loc='best')  # 'loc' can be adjusted to specify the legend position
        # self.spectral_ax.set_xlabel('Wavelength')
        # self.spectral_ax.set_ylabel('Value')
        # self.spectral_ax.set_title('Spectral Plot')

        # xmin, xmax = self.spectral_ax.get_xlim()
        # xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        # xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        # min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
        # buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        # self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

        # self.spectral_canvas = FigureCanvasTkAgg(spectral_figure, master=self.spectral_window)
        # self.spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # # add a button to turn the legend on or off
        # self.toggle_legend_button = tk.Button(ui_frame, text="Legend (On)", command=self.toggle_spectral_plot_legend)
        # self.toggle_legend_button.pack(side=tk.TOP)

        # # Add a button to reset x-axis span
        # self.reset_spectral_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_x_axis_span)
        # self.reset_spectral_x_axis_button.pack(side=tk.TOP)

        # # Add built-in span options to a dropdown menu
        # span_options = ["Full Span", "410 - 1000 nm", "410 - 2500 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
        # self.span_var = tk.StringVar()
        # self.span_var.set("Full Span")  # Set the default span option
        # span_menu = ttk.Combobox(ui_frame, textvariable=self.span_var, values=span_options, state="readonly")
        # span_menu.pack(side=tk.TOP)

        # # Bind an event to update the x-axis span when a span option is selected
        # span_menu.bind("<<ComboboxSelected>>", self.update_x_axis_span)

        # # add a button to ignore bad bands
        # self.ignore_bad_bands_button = tk.Button(ui_frame, text="Ignore Bad Bands (Off)", command=self.toggle_ignore_bad_bands)
        # self.ignore_bad_bands_button.pack(side=tk.TOP)

        # # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        # self.library_spectra_var = tk.StringVar(ui_frame)
        # self.library_spectra_var.set("None") # default value
        # self.library_spectra_list = [name.split('/')[-1] for name in self.usgs_spectra_folders]
        # # sort library_spectra_list alphabetically
        # self.library_spectra_list.sort()
        # # label the librar spectra drop down menu "library spectra"
        # self.library_spectra_label = tk.Label(ui_frame, text="Library Spectra")
        # self.library_spectra_label.pack(side=tk.TOP)
        # # create the drop down menu
        # self.library_spectra_dropdown = tk.OptionMenu(ui_frame, self.library_spectra_var, *self.library_spectra_list, command=self.plot_library_spectra)
        # self.library_spectra_dropdown.pack(side=tk.TOP)

        # # add a button to remove library spectra from the plot
        # self.remove_library_spectra_button = tk.Button(ui_frame, text="Remove Library Spectra", command=self.remove_library_spectra)
        # self.remove_library_spectra_button.pack(side=tk.TOP)

        # # add a button to collect the point and save it to a table
        # self.collect_point_button = tk.Button(ui_frame, text="Collect Point", command=self.collect_point)
        # self.collect_point_button.pack(side=tk.TOP)

        # # Create a toolbar for the spectral plot
        # toolbar = NavigationToolbar2Tk(self.spectral_canvas, self.spectral_window)
        # toolbar.update()
        # self.spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # # Add span selector for x-axis
        # self.create_x_axis_span_selector(self.left_wvl)

    def toggle_spectral_plot_legend(self):
        if self.spectral_ax.get_legend() is None:
            self.spectral_ax.legend()
            self.toggle_legend_button.config(text="Legend (On)", relief="raised")
        else:
            self.spectral_ax.legend_.remove()
            self.toggle_legend_button.config(text="Legend (Off)", relief="sunken")
        self.spectral_canvas.draw()

    def plot_library_spectra(self, event):
        # plot the selected library spectrum
        spec_name = event
        path_to_spec_data = self.usgs_spectra_path + spec_name + "/" + spec_name + ".txt"
        path_to_wvl_data = (
            self.usgs_spectra_path + spec_name + "/" + "*Wavelengths*" + ".txt"
        )
        print(path_to_spec_data)
        library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
        library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
        xmin, xmax = self.spectral_ax.get_xlim()
        ymin, ymax = self.spectral_ax.get_ylim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))

        # Check if spectrum has valid data for scaling
        spectrum_slice = self.spectrum[xmin_idx:xmax_idx]
        if np.all(np.isnan(spectrum_slice)):
            messagebox.showwarning(
                "Warning",
                "Cannot add library spectrum: ROI spectrum contains no valid data."
            )
            return

        # Scale library spectrum to match data range
        scaled_reflectance = (library_reflectance - np.nanmin(library_reflectance)) * (
            (np.nanmax(spectrum_slice) - np.nanmin(spectrum_slice))
            / (np.nanmax(library_reflectance) - np.nanmin(library_reflectance))
        ) + np.nanmin(spectrum_slice)

        # Plot and get line reference
        (line,) = self.spectral_ax.plot(library_wvl, scaled_reflectance, label=spec_name)

        # Store library spectrum data for stretch/offset controls
        # Store the fixed pivot (median) and initialize anchor offset to 0
        pivot_value = np.nanmedian(scaled_reflectance)
        self.library_spectra_data[spec_name] = {
            'line': line,
            'original_y': scaled_reflectance.copy(),
            'wvl': library_wvl,
            'stretch': 1.0,
            'offset': 0.0,
            'anchor': pivot_value,
            'pivot': pivot_value  # Fixed reference point for stretching
        }

        self.spectral_ax.set_ylim((ymin, ymax))
        self.spectral_ax.set_xlim(xmin, xmax)
        self.spectral_ax.legend()
        self.spectral_canvas.draw()

        # Update stretch/offset controls
        self.update_library_spectra_controls()

    def update_library_spectra_controls(self):
        """Create or update the stretch/offset controls for library spectra."""
        # Create controls frame if it doesn't exist
        if self.library_spectra_controls_frame is None or not self.library_spectra_controls_frame.winfo_exists():
            self.library_spectra_controls_frame = tk.Toplevel(self.spectral_window)
            self.library_spectra_controls_frame.title("Library Spectra Controls")
            self.library_spectra_controls_frame.geometry("400x300")

        # Clear existing controls
        for widget in self.library_spectra_controls_frame.winfo_children():
            widget.destroy()

        if not self.library_spectra_data:
            tk.Label(self.library_spectra_controls_frame, text="No library spectra loaded").pack()
            return

        # Create scrollable frame
        canvas = tk.Canvas(self.library_spectra_controls_frame)
        scrollbar = ttk.Scrollbar(self.library_spectra_controls_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Create controls for each library spectrum
        for spec_name, spec_data in self.library_spectra_data.items():
            frame = ttk.LabelFrame(scrollable_frame, text=spec_name, padding=5)
            frame.pack(fill="x", padx=5, pady=5)

            # Stretch controls
            stretch_frame = ttk.Frame(frame)
            stretch_frame.pack(fill="x")
            ttk.Label(stretch_frame, text="Stretch:", width=8).pack(side="left")

            stretch_var = tk.DoubleVar(value=spec_data['stretch'])
            stretch_entry = ttk.Entry(stretch_frame, textvariable=stretch_var, width=8)
            stretch_entry.pack(side="left", padx=2)

            # Expanded range: 0.01 to 10.0 for more flexibility
            stretch_scale = ttk.Scale(stretch_frame, from_=0.01, to=10.0, variable=stretch_var, orient="horizontal", length=120)
            stretch_scale.pack(side="left", padx=2)

            ttk.Button(stretch_frame, text="-", width=2,
                       command=lambda n=spec_name, v=stretch_var: self.adjust_stretch(n, v, -0.05)).pack(side="left")
            ttk.Button(stretch_frame, text="+", width=2,
                       command=lambda n=spec_name, v=stretch_var: self.adjust_stretch(n, v, 0.05)).pack(side="left")

            # Offset controls
            offset_frame = ttk.Frame(frame)
            offset_frame.pack(fill="x")
            ttk.Label(offset_frame, text="Offset:", width=8).pack(side="left")

            offset_var = tk.DoubleVar(value=spec_data['offset'])
            offset_entry = ttk.Entry(offset_frame, textvariable=offset_var, width=8)
            offset_entry.pack(side="left", padx=2)

            # Get y-range for offset scale - use larger range (3x spectrum range)
            y_range = np.nanmax(self.spectrum) - np.nanmin(self.spectrum)
            offset_scale = ttk.Scale(offset_frame, from_=-y_range*3, to=y_range*3, variable=offset_var, orient="horizontal", length=120)
            offset_scale.pack(side="left", padx=2)

            ttk.Button(offset_frame, text="-", width=2,
                       command=lambda n=spec_name, v=offset_var, r=y_range: self.adjust_offset(n, v, -r*0.01)).pack(side="left")
            ttk.Button(offset_frame, text="+", width=2,
                       command=lambda n=spec_name, v=offset_var, r=y_range: self.adjust_offset(n, v, r*0.01)).pack(side="left")

            # Anchor controls
            anchor_frame = ttk.Frame(frame)
            anchor_frame.pack(fill="x")
            ttk.Label(anchor_frame, text="Anchor:", width=8).pack(side="left")

            anchor_var = tk.DoubleVar(value=spec_data['anchor'])
            anchor_entry = ttk.Entry(anchor_frame, textvariable=anchor_var, width=8)
            anchor_entry.pack(side="left", padx=2)

            # Anchor scale based on broader range (2x spectrum range centered on median)
            spec_median = np.nanmedian(spec_data['original_y'])
            spec_range = np.nanmax(spec_data['original_y']) - np.nanmin(spec_data['original_y'])
            anchor_scale = ttk.Scale(anchor_frame, from_=spec_median - spec_range, to=spec_median + spec_range,
                                    variable=anchor_var, orient="horizontal", length=120)
            anchor_scale.pack(side="left", padx=2)

            ttk.Button(anchor_frame, text="-", width=2,
                       command=lambda n=spec_name, v=anchor_var, r=y_range: self.adjust_anchor(n, v, -r*0.01)).pack(side="left")
            ttk.Button(anchor_frame, text="+", width=2,
                       command=lambda n=spec_name, v=anchor_var, r=y_range: self.adjust_anchor(n, v, r*0.01)).pack(side="left")

            # Bind variable changes to update function with error handling
            def safe_update_wrapper(spec_name, stretch_var, offset_var, anchor_var):
                """Wrapper to safely handle variable updates with validation."""
                try:
                    stretch = stretch_var.get()
                    offset = offset_var.get()
                    anchor = anchor_var.get()
                    self.apply_stretch_offset(spec_name, stretch, offset, anchor)
                except (tk.TclError, ValueError):
                    # Ignore errors from incomplete input (e.g., typing "-" or "0.")
                    pass

            stretch_var.trace_add("write", lambda *args, n=spec_name, sv=stretch_var, ov=offset_var, av=anchor_var:
                                  safe_update_wrapper(n, sv, ov, av))
            offset_var.trace_add("write", lambda *args, n=spec_name, sv=stretch_var, ov=offset_var, av=anchor_var:
                                 safe_update_wrapper(n, sv, ov, av))
            anchor_var.trace_add("write", lambda *args, n=spec_name, sv=stretch_var, ov=offset_var, av=anchor_var:
                                 safe_update_wrapper(n, sv, ov, av))

            # Reset button to restore default values
            reset_frame = ttk.Frame(frame)
            reset_frame.pack(fill="x", pady=5)
            ttk.Button(reset_frame, text="Reset to Defaults",
                      command=lambda n=spec_name, sv=stretch_var, ov=offset_var, av=anchor_var, p=spec_data['pivot']:
                      self.reset_library_spectrum(n, sv, ov, av, p)).pack()

    def adjust_stretch(self, spec_name, var, delta):
        """Adjust stretch value by delta."""
        new_val = max(0.01, var.get() + delta)  # Minimum 0.01, no maximum
        var.set(round(new_val, 3))

    def adjust_offset(self, spec_name, var, delta):
        """Adjust offset value by delta."""
        var.set(round(var.get() + delta, 6))

    def adjust_anchor(self, spec_name, var, delta):
        """Adjust anchor value by delta."""
        var.set(round(var.get() + delta, 6))

    def reset_library_spectrum(self, spec_name, stretch_var, offset_var, anchor_var, pivot):
        """Reset library spectrum to default values."""
        stretch_var.set(1.0)
        offset_var.set(0.0)
        anchor_var.set(pivot)

    def apply_stretch_offset(self, spec_name, stretch, offset, anchor):
        """
        Apply stretch and offset to a library spectrum using an anchor point.

        The spectrum is stretched around a fixed pivot point (original median),
        then shifted vertically by the anchor value.
        Formula: new_y = anchor + stretch * (original_y - pivot) + offset

        Parameters:
        -----------
        spec_name : str
            Name of the library spectrum
        stretch : float
            Stretch factor (multiplies deviation from pivot)
        offset : float
            Additional vertical offset applied after stretching
        anchor : float
            Vertical position for the pivot point (increases anchor moves spectrum up)
        """
        if spec_name not in self.library_spectra_data:
            return

        try:
            # Validate inputs
            stretch = float(stretch)
            offset = float(offset)
            anchor = float(anchor)

            spec_data = self.library_spectra_data[spec_name]
            spec_data['stretch'] = stretch
            spec_data['offset'] = offset
            spec_data['anchor'] = anchor

            # Calculate new y data: anchor + stretch * (original - pivot) + offset
            # Pivot is fixed (original median), anchor controls vertical position
            original_y = spec_data['original_y']
            pivot = spec_data['pivot']
            new_y = anchor + stretch * (original_y - pivot) + offset
            spec_data['line'].set_ydata(new_y)

            self.spectral_canvas.draw()
        except (ValueError, TypeError, KeyError) as e:
            # Silently ignore invalid input during typing
            pass

    def plot_mica_library_spectra(self, event):
        # plot the selected MICA library spectrum
        spec_name = "MICA:" + event
        path_to_spec_data = self.mica_spectra_path + event + ".tab"
        print(path_to_spec_data)
        mica_df = pd.read_csv(path_to_spec_data, header=None)
        library_wvl = mica_df[0].values * 1000.0  # Convert to nm
        library_reflectance = mica_df[1].values
        library_reflectance = np.where(
            library_reflectance > 6000, np.nan, library_reflectance
        )

        xmin, xmax = self.spectral_ax.get_xlim()
        ymin, ymax = self.spectral_ax.get_ylim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))

        # Check if spectrum has valid data for scaling
        spectrum_slice = self.spectrum[xmin_idx:xmax_idx]
        if np.all(np.isnan(spectrum_slice)):
            messagebox.showwarning(
                "Warning",
                "Cannot add library spectrum: ROI spectrum contains no valid data."
            )
            return

        # Scale library spectrum to match data range
        scaled_reflectance = (library_reflectance - np.nanmin(library_reflectance)) * (
            (np.nanmax(spectrum_slice) - np.nanmin(spectrum_slice))
            / (np.nanmax(library_reflectance) - np.nanmin(library_reflectance))
        ) + np.nanmin(spectrum_slice)

        # Plot and get line reference
        (line,) = self.spectral_ax.plot(library_wvl, scaled_reflectance, label=event)

        # Store library spectrum data for stretch/offset controls
        # Store the fixed pivot (median) and initialize anchor offset to 0
        pivot_value = np.nanmedian(scaled_reflectance)
        self.library_spectra_data[spec_name] = {
            'line': line,
            'original_y': scaled_reflectance.copy(),
            'wvl': library_wvl,
            'stretch': 1.0,
            'offset': 0.0,
            'anchor': pivot_value,
            'pivot': pivot_value  # Fixed reference point for stretching
        }

        self.spectral_ax.set_ylim((ymin, ymax))
        self.spectral_ax.set_xlim(xmin, xmax)
        self.spectral_ax.legend()
        self.spectral_canvas.draw()

        # Update stretch/offset controls
        self.update_library_spectra_controls()

    def remove_library_spectra(self):
        # Remove all lines except the first one (assuming it's the main line of the plot)
        for line in self.spectral_ax.lines[1:]:
            line.remove()

        # Clear library spectra tracking data
        self.library_spectra_data = {}

        # Close controls window if open
        if self.library_spectra_controls_frame is not None and self.library_spectra_controls_frame.winfo_exists():
            self.library_spectra_controls_frame.destroy()

        # Update the plot
        self.update_spectral_plot()

    def toggle_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_bad_bands_flag:
                self.ignore_bad_bands_flag = False
                self.ignore_bad_bands_button.config(
                    text="Ignore Bad Bands (Off)", relief="sunken"
                )
            else:
                self.ignore_bad_bands_flag = True
                self.ignore_bad_bands_button.config(
                    text="Ignore Bad Bands (On)", relief="raised"
                )
            self.update_spectral_plot()
        else:
            messagebox.showwarning(
                "Warning",
                "No bad bands set. Set bad bands first. (Processing > Select Bad Bands)",
            )

    def toggle_1um_offset_correction(self):
        """Toggle the 1µm offset correction for CRISM spectra."""
        if self.offset_correction_applied:
            # Turn off correction
            self.offset_correction_applied = False
            self.correct_offset_button.config(
                text="Correct 1µm Offset (Off)", relief="sunken"
            )
        else:
            # Turn on correction
            self.offset_correction_applied = True
            self.correct_offset_button.config(
                text="Correct 1µm Offset (On)", relief="raised"
            )
        self.update_spectral_plot()

    def apply_1um_offset_correction(self, spectrum):
        """
        Correct the 1µm offset in CRISM spectra.

        This corrects the discontinuity between short-wave and long-wave detectors
        around 1000nm by adjusting the short-wavelength values to align with the
        long-wavelength values.

        Parameters:
        -----------
        spectrum : np.ndarray
            The input spectrum to correct

        Returns:
        --------
        corrected_spectrum : np.ndarray
            The corrected spectrum
        """
        wvl = np.array(self.left_wvl)

        # Find the index closest to 1000nm (1µm)
        detector_boundary_idx = np.argmin(np.abs(wvl - 1000))

        # Define windows around the detector boundary
        # Use 5 samples before and after the boundary
        window_size = 5

        # Get shortwave side (just before the boundary)
        sw_start = max(0, detector_boundary_idx - window_size)
        sw_end = detector_boundary_idx

        # Get longwave side (just after the boundary)
        lw_start = detector_boundary_idx + 1
        lw_end = min(len(spectrum), detector_boundary_idx + window_size + 1)

        # Calculate medians on each side of the boundary
        sw_values = spectrum[sw_start:sw_end]
        lw_values = spectrum[lw_start:lw_end]

        # Only calculate if we have valid data on both sides
        if len(sw_values) > 0 and len(lw_values) > 0:
            # Use nanmedian to ignore NaN values
            sw_median = np.nanmedian(sw_values)
            lw_median = np.nanmedian(lw_values)

            # Check if we got valid medians
            if not np.isnan(sw_median) and not np.isnan(lw_median):
                # Calculate the offset
                offset = lw_median - sw_median

                # Apply correction: add offset to all shortwave values
                corrected_spectrum = spectrum.copy()
                corrected_spectrum[:detector_boundary_idx] += offset

                return corrected_spectrum

        # If correction couldn't be applied, return original spectrum
        return spectrum

    def update_spectral_plot(self):
        """
        Update the spectral plot with the current spectrum.

        Uses the new plotting system: sends spectrum to the active plot,
        or creates a new plot if none exists.
        """
        # Clean up any closed windows first
        self.plot_manager.cleanup_closed_windows()

        # Store original spectrum for any correction toggles
        self.original_spectrum = self.spectrum.copy()

        # Check if we have an active plot
        active_plot = self.plot_manager.get_active_plot()

        if active_plot is None:
            # No active plot - create a new one
            print("creating spectral plot")
            self.create_spectral_plot()
        else:
            # Send spectrum to the active plot
            spectrum_data = SpectrumData(
                wavelengths=np.array(self.left_wvl),
                values=self.spectrum.copy(),
                label=self.spectrum_label,
            )

            # Check accumulate mode - if enabled, add without clearing
            # Otherwise, clear previous spectra (maintaining old behavior)
            if not active_plot.is_accumulate_mode():
                active_plot.clear_spectra(update_view=False)
            active_plot.add_spectrum(spectrum_data)

            # Update backwards-compatibility references
            self.spectral_window = active_plot.window
            self.spectral_ax = active_plot._ax
            self.spectral_canvas = active_plot._canvas
            if active_plot._spectrum_lines:
                self.spectral_line = active_plot._spectrum_lines[0]
            self._current_point_plot_window = active_plot

    def _add_roi_spectra_to_active_plot(self) -> None:
        """
        Add ROI spectra to the active spectral plot window.

        This provides a unified approach to plotting both point and ROI spectra
        through the same system.
        """
        # Clean up any closed windows first
        self.plot_manager.cleanup_closed_windows()

        # Get or create active plot
        active_plot = self.plot_manager.get_active_plot()
        if active_plot is None:
            active_plot = self._create_spectral_plot_window(title="ROI Spectral Plot")

        # Check accumulate mode
        if not active_plot.is_accumulate_mode():
            active_plot.clear_spectra(update_view=False)

        # Add each ROI spectrum
        for i, (color, spectrum) in enumerate(zip(self.polygon_colors, self.polygon_spectra)):
            # Get ROI name from table data if available
            if i < len(self.polygons_table_data):
                label = self.polygons_table_data[i][0]  # First column is ROI name
            else:
                label = f"ROI {i+1}"

            spectrum_data = SpectrumData(
                wavelengths=np.array(self.left_wvl),
                values=np.array(spectrum).flatten(),
                label=label,
                color=color,
                metadata={'source': 'roi', 'roi_index': i}
            )
            active_plot.add_spectrum(spectrum_data, update_view=False)

        active_plot._update_view()

    def reset_x_axis_span(self):
        # Reset x-axis span to the default range
        self.spectral_ax.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))

        # Check if the spectrum has any valid data in the displayed range
        spectrum_slice = self.spectrum[xmin_idx:xmax_idx]
        if not np.all(np.isnan(spectrum_slice)):
            min_y, max_y = np.nanmin(spectrum_slice), np.nanmax(spectrum_slice)
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
        self.spectral_canvas.draw()

    def update_x_axis_span(self, event):
        selected_span = self.span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_wvl[0], self.left_wvl[-1]),
            "410 - 1000 nm": (410, 1000),
            "410 - 2500 nm": (410, 2500),
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900),
        }
        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            self.spectral_ax.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))

            # Calculate y-limits based on the data within the new span
            y_min = np.nanmin(self.spectrum[xmin_idx:xmax_idx])
            y_max = np.nanmax(self.spectrum[xmin_idx:xmax_idx])

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            self.spectral_ax.set_ylim(y_min - buffer, y_max + buffer)

            self.spectral_canvas.draw()

    def create_x_axis_span_selector(self, x_data):
        def on_x_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(x_data - xmin))
            xmax_idx = np.argmin(np.abs(x_data - xmax))
            self.spectral_ax.set_xlim(x_data[xmin_idx], x_data[xmax_idx])

            # Calculate y-limits based on the data within the new span
            y_min = np.nanmin(self.spectrum[xmin_idx:xmax_idx])
            y_max = np.nanmax(self.spectrum[xmin_idx:xmax_idx])

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            self.spectral_ax.set_ylim(y_min - buffer, y_max + buffer)

            self.spectral_canvas.draw()

        self.x_span_selector = SpanSelector(
            self.spectral_ax, on_x_span_select, "horizontal", useblit=True
        )

    # ----------------------------------------------------------------
    # ratio point plots
    # ----------------------------------------------------------------
    def create_ratio_spectral_plot(self):
        self.ratio_spectral_window = tk.Toplevel(self.root)
        self.ratio_spectral_window.title("Ratio Spectral Plot")

        # Create a frame to hold UI elements with a fixed size
        ratio_ui_frame = tk.Frame(self.ratio_spectral_window)
        ratio_ui_frame.pack(fill=tk.X)

        # ratio_spectral_figure, self.ratio_spectral_ax = plt.subplots(figsize=(5,3))
        self.ratio_spectral_figure, (
            self.ratio_spectral_ax1,
            self.ratio_spectral_ax2,
        ) = plt.subplots(2, 1, figsize=(5, 6))
        self.ratio_spectral_axes = (self.ratio_spectral_ax1, self.ratio_spectral_ax2)

        # numerator and denominator are plotted on the top
        (self.numerator_spectral_line,) = self.ratio_spectral_ax1.plot(
            self.left_wvl, self.spectrum, label="numerator"
        )
        (self.denominator_spectral_line,) = self.ratio_spectral_ax1.plot(
            self.left_wvl, self.denom_spectrum, label="denominator"
        )
        xmin, xmax = self.ratio_spectral_ax1.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(
            (
                np.nanmin(self.spectrum[xmin_idx:xmax_idx]),
                np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]),
            )
        ), np.nanmax(
            (
                np.nanmax(self.spectrum[xmin_idx:xmax_idx]),
                np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx]),
            )
        )
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.ratio_spectral_ax1.set_ylim(min_y - buffer, max_y + buffer)

        # ratioed spectrum plots on the bottom
        (self.ratio_spectral_line,) = self.ratio_spectral_ax2.plot(
            self.left_wvl, self.ratio_spectrum, label="ratio"
        )
        xmin, xmax = self.ratio_spectral_ax2.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(
            self.ratio_spectrum[xmin_idx:xmax_idx]
        )
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.ratio_spectral_ax2.set_ylim(min_y - buffer, max_y + buffer)

        self.ratio_spectral_ax1.legend(
            loc="best"
        )  # 'loc' can be adjusted to specify the legend position
        self.ratio_spectral_ax2.legend(
            loc="best"
        )  # 'loc' can be adjusted to specify the legend position

        self.ratio_spectral_canvas = FigureCanvasTkAgg(
            self.ratio_spectral_figure, master=self.ratio_spectral_window
        )
        self.ratio_spectral_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # ---
        # add a button to turn the legend on or off
        self.toggle_ratio_legend_button = tk.Button(
            ratio_ui_frame,
            text="Legend (On)",
            command=self.toggle_ratio_spectral_plot_legend,
        )
        self.toggle_ratio_legend_button.pack(side=tk.TOP)

        # add a button to ignore bad bands
        self.ignore_bad_bands_button = tk.Button(
            ratio_ui_frame,
            text="Ignore Bad Bands (Off)",
            command=self.toggle_ratio_ignore_bad_bands,
        )
        self.ignore_bad_bands_button.pack(side=tk.TOP)

        # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        self.ratio_library_spectra_var = tk.StringVar(ratio_ui_frame)
        self.ratio_library_spectra_var.set("None")  # default value
        self.ratio_library_spectra_list = [
            name.split("/")[-1] for name in self.usgs_spectra_folders
        ]
        # sort library_spectra_list alphabetically
        self.ratio_library_spectra_list.sort()
        # label the library spectra drop down menu "library spectra"
        self.ratio_library_spectra_label = tk.Label(
            ratio_ui_frame, text="Library Spectra"
        )
        self.ratio_library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.ratio_library_spectra_dropdown = tk.OptionMenu(
            ratio_ui_frame,
            self.ratio_library_spectra_var,
            *self.ratio_library_spectra_list,
            command=self.plot_ratio_library_spectra,
        )
        self.ratio_library_spectra_dropdown.pack(side=tk.TOP)

        # add a button to remove library spectra from the plot
        self.remove_ratio_library_spectra_button = tk.Button(
            ratio_ui_frame,
            text="Remove Library Spectra",
            command=self.remove_ratio_library_spectra,
        )
        self.remove_ratio_library_spectra_button.pack(side=tk.TOP)

        # add a button to collect the point and save it to a table
        self.collect_ratio_point_button = tk.Button(
            ratio_ui_frame, text="Collect Point", command=self.collect_point
        )
        self.collect_ratio_point_button.pack(side=tk.TOP)
        # ----

        # Add a button to reset x-axis span
        self.reset_ratio_x_axis_button = tk.Button(
            ratio_ui_frame,
            text="Reset X-Axis Span",
            command=self.reset_ratio_x_axis_span,
        )
        self.reset_ratio_x_axis_button.pack(side=tk.RIGHT)

        # Add built-in span options to a dropdown menu
        span_options = [
            "Full Span",
            "410 - 1000 nm",
            "410 - 2500 nm",
            "1000 - 2600 nm",
            "1200 - 2000 nm",
            "1800 - 2500 nm",
            "2000 - 2500 nm",
            "2700 - 3900 nm",
        ]
        self.ratio_span_var = tk.StringVar()
        self.ratio_span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(
            ratio_ui_frame,
            textvariable=self.ratio_span_var,
            values=span_options,
            state="readonly",
        )
        span_menu.pack(side=tk.RIGHT)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_ratio_x_axis_span)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(
            self.ratio_spectral_canvas, self.ratio_spectral_window
        )
        toolbar.update()
        self.ratio_spectral_canvas.get_tk_widget().pack(
            side=tk.TOP, fill=tk.BOTH, expand=True
        )

        # Add span selector for x-axis
        self.create_ratio_x_axis_span_selector(self.left_wvl)

    def update_ratio_spectral_plot(self):
        if (
            self.ratio_spectral_window is None
            or not self.ratio_spectral_window.winfo_exists()
        ):
            self.create_ratio_spectral_plot()
        else:

            # numerator and denominator are plotted on the top
            self.numerator_spectral_line.set_ydata(self.spectrum)
            self.denominator_spectral_line.set_ydata(self.denom_spectrum)
            xmin, xmax = self.ratio_spectral_ax1.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin(
                (
                    np.nanmin(self.spectrum[xmin_idx:xmax_idx]),
                    np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]),
                )
            ), np.nanmax(
                (
                    np.nanmax(self.spectrum[xmin_idx:xmax_idx]),
                    np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx]),
                )
            )
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax1.set_ylim(min_y - buffer, max_y + buffer)

            # ratioed spectrum plots on the bottom
            self.ratio_spectral_line.set_ydata(self.ratio_spectrum)
            xmin, xmax = self.ratio_spectral_ax2.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(
                self.ratio_spectrum[xmin_idx:xmax_idx]
            )
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax2.set_ylim(min_y - buffer, max_y + buffer)

            self.ratio_spectral_ax1.legend(
                loc="best"
            )  # 'loc' can be adjusted to specify the legend position
            self.ratio_spectral_ax2.legend(
                loc="best"
            )  # 'loc' can be adjusted to specify the legend position

            self.ratio_spectral_canvas.draw()

    def plot_ratio_library_spectra(self, event):
        # plot the selected library spectrum
        # self.point_plot_library_indices
        path_to_spec_data = self.usgs_spectra_path + event + "/" + event + ".txt"
        path_to_wvl_data = (
            self.usgs_spectra_path + event + "/" + "*Wavelengths*" + ".txt"
        )
        print(path_to_spec_data)
        library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
        library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
        for i, ax in enumerate(self.ratio_spectral_axes):
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            # plot the values in spectrum_df on ax
            ax.plot(
                library_wvl,
                library_reflectance,
                label=self.ratio_library_spectra_var.get(),
            )
            # scale each line so that the min and max values equal the min and max max value of self.spectrum
            for line in ax.lines[1:]:
                if i == 0:
                    s = self.spectrum
                elif i == 1:
                    s = self.ratio_spectrum
                new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * (
                    (np.nanmax(s[xmin_idx:xmax_idx]) - np.nanmin(s[xmin_idx:xmax_idx]))
                    / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))
                ) + np.nanmin(s[xmin_idx:xmax_idx])
                line.set_ydata(new_y_data)
            # min_y, max_y = np.nanmin(library_reflectance), np.nanmax(library_reflectance)
            # buffer = (max_y - min_y) * 0.1
            # new_y_lim = (np.nanmin((ymin, min_y - buffer)), np.nanmax((ymax, max_y + buffer)))
            ax.set_ylim((ymin, ymax))
            ax.set_xlim(xmin, xmax)
            ax.legend()
            self.ratio_spectral_canvas.draw()

    def remove_ratio_library_spectra(self):
        # Remove all lines except the first one (assuming it's the main line of the plot)
        for line in self.ratio_spectral_ax1.lines[2:]:
            line.remove()
        for line in self.ratio_spectral_ax2.lines[1:]:
            line.remove()

        # Update the plot
        self.update_ratio_spectral_plot()

    def toggle_ratio_spectral_plot_legend(self):
        for ax in self.ratio_spectral_axes:
            if ax.get_legend() is None:
                ax.legend()
                self.toggle_ratio_legend_button.config(
                    text="Legend (On)", relief="raised"
                )
            else:
                ax.legend_.remove()
                self.toggle_ratio_legend_button.config(
                    text="Legend (Off)", relief="sunken"
                )
            self.ratio_spectral_canvas.draw()

    def toggle_ratio_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_ratio_bad_bands_flag:
                self.ignore_ratio_bad_bands_flag = False
                self.ignore_ratio_bad_bands_flag.config(
                    text="Ignore Bad Bands (Off)", relief="sunken"
                )
            else:
                self.ignore_ratio_bad_bands_flag = True
                self.ignore_ratio_bad_bands_flag.config(
                    text="Ignore Bad Bands (On)", relief="raised"
                )
            self.update_spectral_plot()
        else:
            messagebox.showwarning(
                "Warning",
                "No bad bands set. Set bad bands first. (Processing > Select Bad Bands)",
            )

    def reset_ratio_x_axis_span(self):
        # Reset x-axis span to the default range
        self.ratio_spectral_ax1.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.ratio_spectral_ax1.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(
            (
                np.nanmin(self.spectrum[xmin_idx:xmax_idx]),
                np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]),
            )
        ), np.nanmax(
            (
                np.nanmax(self.spectrum[xmin_idx:xmax_idx]),
                np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx]),
            )
        )
        buffer = (max_y - min_y) * 0.1
        self.ratio_spectral_ax1.set_ylim(min_y - buffer, max_y + buffer)

        self.ratio_spectral_ax2.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.ratio_spectral_ax2.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(
            self.ratio_spectrum[xmin_idx:xmax_idx]
        )
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.ratio_spectral_ax2.set_ylim(min_y - buffer, max_y + buffer)

        self.ratio_spectral_canvas.draw()

    def update_ratio_x_axis_span(self, event):
        selected_span = self.ratio_span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_wvl[0], self.left_wvl[-1]),
            "410 - 1000 nm": (410, 1000),
            "410 - 2500 nm": (410, 2500),
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900),
        }
        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            self.ratio_spectral_ax1.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))
            min_y, max_y = np.nanmin(
                (
                    np.nanmin(self.spectrum[xmin_idx:xmax_idx]),
                    np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]),
                )
            ), np.nanmax(
                (
                    np.nanmax(self.spectrum[xmin_idx:xmax_idx]),
                    np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx]),
                )
            )
            buffer = (max_y - min_y) * 0.1
            self.ratio_spectral_ax1.set_ylim(min_y - buffer, max_y + buffer)

            self.ratio_spectral_ax2.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))

            # Calculate y-limits based on the data within the new span
            y_min = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx])
            y_max = np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax2.set_ylim(y_min - buffer, y_max + buffer)

            self.ratio_spectral_canvas.draw()

    def create_ratio_x_axis_span_selector(self, x_data):
        def on_x_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(x_data - xmin))
            xmax_idx = np.argmin(np.abs(x_data - xmax))
            self.ratio_spectral_ax1.set_xlim(x_data[xmin_idx], x_data[xmax_idx])
            self.ratio_spectral_ax2.set_xlim(x_data[xmin_idx], x_data[xmax_idx])

            # Calculate y-limits based on the data within the new span
            y_min, y_max = np.nanmin(
                (
                    np.nanmin(self.spectrum[xmin_idx:xmax_idx]),
                    np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]),
                )
            ), np.nanmax(
                (
                    np.nanmax(self.spectrum[xmin_idx:xmax_idx]),
                    np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx]),
                )
            )
            buffer = (y_max - y_min) * 0.1
            self.ratio_spectral_ax1.set_ylim(y_min - buffer, y_max + buffer)

            y_min = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx])
            y_max = np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax2.set_ylim(y_min - buffer, y_max + buffer)

            self.ratio_spectral_canvas.draw()

        self.ratio_top_x_span_selector = SpanSelector(
            self.ratio_spectral_ax1, on_x_span_select, "horizontal", useblit=True
        )

        self.ratio_btm_x_span_selector = SpanSelector(
            self.ratio_spectral_ax2, on_x_span_select, "horizontal", useblit=True
        )

    # ----------------------------------------------------------------
    # Spectral Processing
    # ----------------------------------------------------------------
    def select_bad_bands(self):
        # function to interactively select bad bands
        # create a figure with the spectral plot
        self.bad_bands_window = tk.Toplevel(self.root)
        self.bad_bands_window.title("Bad Bands Selector")

        # grab a spectrum
        s1, s2, s3 = np.shape(self.left_data)
        self.bbl_spectrum = self.left_data[
            int(
                s1 / 2,
            ),
            int(
                s2 / 2,
            ),
            :,
        ].flatten()

        # create the spectral plot
        bad_bands_figure, self.bad_bands_ax = plt.subplots(figsize=(5, 3))
        (self.bad_bands_line,) = self.bad_bands_ax.plot(
            self.left_wvl, self.bbl_spectrum
        )
        self.bad_bands_ax.set_xlabel("Wavelength")
        self.bad_bands_ax.set_ylabel("Value")
        self.bad_bands_ax.set_title("Spectral Plot")
        xmin, xmax = self.bad_bands_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.bbl_spectrum[xmin_idx:xmax_idx]), np.nanmax(
            self.bbl_spectrum[xmin_idx:xmax_idx]
        )
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.bad_bands_ax.set_ylim(min_y - buffer, max_y + buffer)

        # place the plot in the bad_bands_window
        self.bad_bands_canvas = FigureCanvasTkAgg(
            bad_bands_figure, master=self.bad_bands_window
        )
        self.bad_bands_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # initialize the bad_bands list as all 1's (i.e., good)
        if hasattr(self, "bad_bands"):
            if self.bad_bands is not None:
                self.bad_bands_line.set_ydata(
                    np.where(np.array(self.bad_bands) == 0, np.nan, self.bbl_spectrum)
                )
                self.bad_bands_canvas.draw()
            else:
                self.bad_bands = [1] * len(self.left_wvl)
        else:
            self.bad_bands = [1] * len(self.left_wvl)

        # add a span selector to retrieve the range of bad bands
        def on_x_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(self.left_wvl - xmin))
            xmax_idx = np.argmin(np.abs(self.left_wvl - xmax))
            # set the values within the selected range to zero
            self.bad_bands[xmin_idx:xmax_idx] = [0] * (xmax_idx - xmin_idx)
            # update the plot so that the bad bands are nan
            self.bad_bands_line.set_ydata(
                np.where(np.array(self.bad_bands) == 0, np.nan, self.bbl_spectrum)
            )
            self.bad_bands_canvas.draw()
            # update the bad band values in the table
            for i, value in enumerate(self.bad_bands):
                self.bad_bands_table.set(i, column="#2", value=value)

        self.bad_bands_x_span_selector = SpanSelector(
            self.bad_bands_ax, on_x_span_select, "horizontal", useblit=True
        )

        # add a way to show the bad bands in a table
        self.bad_bands_table = ttk.Treeview(
            self.bad_bands_window, columns=("wavelength", "bad")
        )
        self.bad_bands_table.heading("#0", text="Band")
        self.bad_bands_table.heading("wavelength", text="Wavelength")
        self.bad_bands_table.heading("bad", text="Bad=0")
        self.bad_bands_table.pack(fill=tk.BOTH, expand=True)

        # populate the table with the bad bands
        for i, value in enumerate(self.bad_bands):
            self.bad_bands_table.insert(
                parent="",
                index="end",
                iid=i,
                text=str(i + 1),
                values=(self.left_wvl[i], value),
            )

        # add a way to change the bad bands
        def on_bad_band_change(event):
            selection = self.bad_bands_table.selection()
            if selection:  # Check if there is a selection
                item = selection[0]
                column = self.bad_bands_table.identify_column(event.x)
                if column == "#2":  # bad column
                    wvl, v = self.bad_bands_table.item(item)["values"]
                    if v == 1:
                        self.bad_bands_table.set(item, column="#2", value=0)
                    elif v == 0:
                        self.bad_bands_table.set(item, column="#2", value=1)
            # update self.bad_bands to reflect changes in the table
            self.bad_bands = [
                int(self.bad_bands_table.item(i)["values"][1])
                for i in range(len(self.bad_bands))
            ]
            # update the spectral plot to reflect the change in bad bands list
            self.bad_bands_line.set_ydata(
                np.where(np.array(self.bad_bands) == 0, np.nan, self.bbl_spectrum)
            )
            self.bad_bands_canvas.draw()

        self.bad_bands_table.bind("<Double-1>", on_bad_band_change)

        # add an option to reset the bad bands list to all good bands
        self.plot_spectrum_from_main_view = tk.Button(
            self.bad_bands_window,
            text="Plot Spectrum from Main View",
            command=self.get_main_view_spec,
        )
        self.plot_spectrum_from_main_view.pack(side=tk.BOTTOM)

        # add an option to reset the bad bands list to all good bands
        self.reset_bad_bands_button = tk.Button(
            self.bad_bands_window, text="Reset Bad Bands", command=self.reset_bad_bands
        )
        self.reset_bad_bands_button.pack(side=tk.BOTTOM)

        # add an option to linearly interpolate across the bad bands for the whole image cube
        self.interpolate_bad_bands_button = tk.Button(
            self.bad_bands_window,
            text="Interpolate Bad Bands",
            command=self.interpolate_bad_bands,
        )
        self.interpolate_bad_bands_button.pack(side=tk.BOTTOM)

        # add an option to restore the original cube
        self.restore_original_cube = tk.Button(
            self.bad_bands_window,
            text="Restore Uninterpolated Cube",
            command=self.restore_original_cube,
        )
        self.restore_original_cube.pack(side=tk.BOTTOM)

    def reset_bad_bands(self):
        self.bad_bands = [1] * len(self.left_wvl)
        # update the spectral plot to reflect the change in bad bands list
        self.bad_bands_line.set_ydata(
            np.where(np.array(self.bad_bands) == 0, np.nan, self.bbl_spectrum)
        )
        self.bad_bands_canvas.draw()
        # update the bad band values in the table
        for i, value in enumerate(self.bad_bands):
            self.bad_bands_table.set(i, column="#2", value=value)

    def get_main_view_spec(self):
        # get the spectrum from the main view and plot it in the bad_bands_window
        if not self.spectrum.tolist():
            messagebox.showwarning(
                "Warning",
                "No spectrum selected, click anywhere on the main view to load a spectrum.",
            )
        else:
            self.bbl_spectrum = self.spectrum
            self.bad_bands_line.set_ydata(self.spectrum)
            xmin, xmax = self.bad_bands_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin(self.bbl_spectrum[xmin_idx:xmax_idx]), np.nanmax(
                self.bbl_spectrum[xmin_idx:xmax_idx]
            )
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.bad_bands_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.bad_bands_canvas.draw()

    def interpolate_bad_bands(self):
        self.left_data_og = (
            self.left_data
        )  # store the original cube in case we need it.
        interp_cube = np.array(self.ld)
        # set bad bands to nan
        for i, band in enumerate(self.bad_bands):
            if band == 0:
                interp_cube[:, :, i] = np.nan * interp_cube[:, :, i]
            elif band == 1:
                pass

        print(np.shape(interp_cube))

        ni = interpNaNs(interp_cube, self.left_wvl)
        ni.linearInterp()
        self.left_data = ni.data_cube

    def restore_original_cube(self):
        if hasattr(self, "left_data_og"):
            self.left_data = self.left_data_og
        else:
            messagebox.showwarning(
                "Warning", "No original cube found. Run interpolation first."
            )

    def calculate_spectral_parameters(self):
        self.spectral_parameters_window = tk.Toplevel(self.root)
        self.spectral_parameters_window.title("Spectral Parameterization")

        # Checkbox variables
        self.crop_var = tk.IntVar()
        self.bbl_var = tk.IntVar()
        self.interpNans_var = tk.IntVar()
        self.denoise_var = tk.IntVar()

        # Checkboxes
        tk.Checkbutton(
            self.spectral_parameters_window, text="Crop", variable=self.crop_var
        ).pack()
        # Tooltip(self.spectral_parameters_window, "Crop: Specify crop region, like [row0, row1, column0, column1], with row and column values in pixel coordinates (not lat/lon)")

        tk.Checkbutton(
            self.spectral_parameters_window,
            text="BBL",
            variable=self.bbl_var,
            command=self.activate_bbl_function,
        ).pack()
        # Tooltip(self.spectral_parameters_window, "BBL: Bad Bands List, 1=good, 0=bad. Use Processing>Select Bad Bands to select bad bands.")

        tk.Checkbutton(
            self.spectral_parameters_window,
            text="Interpolate NaNs",
            variable=self.interpNans_var,
        ).pack()
        # Tooltip(self.spectral_parameters_window, "Interpolate NaNs: Option to interpolate NaNs")

        tk.Checkbutton(
            self.spectral_parameters_window, text="Denoise", variable=self.denoise_var
        ).pack()
        # Tooltip(self.spectral_parameters_window, "Denois: Option to denoise, this may take a long time. Recommend using with Interpolate NaNs")

        # run button
        tk.Button(
            self.spectral_parameters_window,
            text="Run",
            command=self.run_spectral_parameters,
        ).pack()
        # run button
        tk.Button(
            self.spectral_parameters_window,
            text="Initiate",
            command=self.initiate_spectral_parameters,
        ).pack()

    def initiate_spectral_parameters(self):
        # Instantiate the class and select your input image and output directory
        print(self.bbl_var)
        if self.bbl_var.get() == 1:
            bbl = self.bad_bands
        else:
            bbl = [None]
        if self.interpNans_var.get() == 1:
            interpNans = True
        else:
            interpNans = False
        if self.crop_var.get() == 1:
            crop = self.get_crop_limits()
        else:
            crop = None
        if self.denoise_var.get() == 1:
            denoise = True
        else:
            denoise = False
        self.pc = None
        self.pc = cubeParamCalculator(
            bbl=bbl, interpNans=interpNans, crop=crop, denoise=denoise
        )
        self.edit_valid_params()

    def run_spectral_parameters(self):
        # Run the calculator and save the results
        self.pc.run()

    def get_crop_limits(self):
        # have user input top, bottom, left, right limits for the crop parameters
        # pop up a window with 4 entry boxes properly labeled
        crop_window = Toplevel(self.root)
        crop_window.title("Crop Parameters")
        # Create labels and entry boxes for crop limits
        top_label = tk.Label(crop_window, text="Top:")
        top_label.grid(row=0, column=0, padx=5, pady=5)
        top_entry = Entry(crop_window)
        top_entry.grid(row=0, column=1, padx=5, pady=5)

        bottom_label = tk.Label(crop_window, text="Bottom:")
        bottom_label.grid(row=1, column=0, padx=5, pady=5)
        bottom_entry = Entry(crop_window)
        bottom_entry.grid(row=1, column=1, padx=5, pady=5)

        left_label = tk.Label(crop_window, text="Left:")
        left_label.grid(row=2, column=0, padx=5, pady=5)
        left_entry = Entry(crop_window)
        left_entry.grid(row=2, column=1, padx=5, pady=5)

        right_label = tk.Label(crop_window, text="Right:")
        right_label.grid(row=3, column=0, padx=5, pady=5)
        right_entry = Entry(crop_window)
        right_entry.grid(row=3, column=1, padx=5, pady=5)

        # Function to retrieve the values from the entry boxes
        def get_crop_limits_():
            top = int(top_entry.get())
            bottom = int(bottom_entry.get())
            left = int(left_entry.get())
            right = int(right_entry.get())
            crop_window.destroy()
            return [top, bottom, left, right]

        # Create a button to submit the crop limits
        submit_button = Button(crop_window, text="Submit", command=get_crop_limits_)
        submit_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

        # return the values in a list
        return get_crop_limits_()

    def edit_valid_params(self):
        # Open a new window for editing valid_params
        edit_window = Toplevel()
        edit_window.title("Edit Valid Parameters")

        # Display current valid_params list
        valid_params_listbox = Listbox(edit_window)
        for param in self.pc.validParams:
            valid_params_listbox.insert(END, param)
        valid_params_listbox.pack()

        # Function to add a parameter to valid_params
        def add_parameter():
            param = input_entry.get()
            self.pc.validParams.append(param)
            valid_params_listbox.insert(END, param)

        # Entry widget to input new parameter
        input_entry = Entry(edit_window)
        input_entry.pack()

        # Button to add parameter
        add_button = Button(edit_window, text="Add Parameter", command=add_parameter)
        add_button.pack()

        # Function to remove selected parameter from valid_params
        def remove_parameter():
            index = valid_params_listbox.curselection()[0]
            param = valid_params_listbox.get(index)
            self.pc.validParams.remove(param)
            valid_params_listbox.delete(index)

        # Button to remove selected parameter
        remove_button = Button(
            edit_window, text="Remove Parameter", command=remove_parameter
        )
        remove_button.pack()

        # Button to close the edit window
        close_button = Button(edit_window, text="Close", command=edit_window.destroy)
        close_button.pack()

    def activate_bbl_function(self):
        if not hasattr(self, "bad_bands"):
            self.select_bad_bands()

    def apply_mnf(self):
        # apply MNF to the left data
        data = self.ld
        signal = spectral.calc_stats(data)
        noise = spectral.noise_from_diffs(data)
        self.mnf = spectral.mnf(signal, noise)

    # ----------------------------------------------------------------
    # saving and closing functions
    # ----------------------------------------------------------------
    def save_state(self):
        # Ask the user to choose the file name and location
        file_path = filedialog.asksaveasfilename(
            defaultextension=".pkl",  # Default extension for your serialized data
            filetypes=[("Pickle files", "*.pkl")],  # Filter for pickle files
        )

        if file_path:
            # Create a dictionary to store the relevant attributes
            state_dict = {
                "left_file": self.left_data.filename.replace(".img", ".hdr"),
                "right_file": self.right_data.filename.replace(".img", ".hdr"),
                "all_polygons": self.all_polygons,
                "polygon_colors": self.polygon_colors,
                "polygon_spectra": self.polygon_spectra,
                "stored_ratio_spectra": self.ratio_spectrum_list,
                "best_denom_ids": self.best_denom_id_list,
                "polygon_params": self.polygon_params,
                "polyong_table_data": self.polygons_table_data,
                "spectrum": self.spectrum,
                "ratio_spectrum": self.ratio_spectrum,
                "bad_bands": self.bad_bands,
            }

            # Serialize and save the state dictionary to the chosen file
            with open(file_path, "wb") as file:
                pickle.dump(state_dict, file)

    def load_state(self):
        # Ask the user to choose the saved session file
        file_path = filedialog.askopenfilename(
            filetypes=[("Pickle files", "*.pkl")],  # Filter for pickle files
        )

        if file_path:
            # Load the saved state from the chosen file
            with open(file_path, "rb") as file:
                restored_instance = pickle.load(file)

            # Update the current instance with the restored state
            self.load_left_data(restored_instance["left_file"])
            self.load_right_data(restored_instance["right_file"])
            self.all_polygons = restored_instance["all_polygons"]
            self.polygon_colors = restored_instance["polygon_colors"]
            self.polygon_spectra = restored_instance["polygon_spectra"]
            self.polygon_params = restored_instance["polygon_params"]
            self.spectrum = restored_instance["spectrum"]
            self.polygons_table_data = restored_instance["polyong_table_data"]
            self.ratio_spectrum = restored_instance["ratio_spectrum"]
            self.bad_bands = restored_instance["bad_bands"]
            if self.all_polygons:
                self.create_polygons_menu_window()
                self.create_polygons_table()
                self.update_polygons_table()

    def on_closing(self):
        self.root.destroy()
        self.root.quit


# ----------------------------------------------------------------
# Spectral window class
# ----------------------------------------------------------------
class spectral_window(SpectralCubeAnalysisTool):
    def __init__(self, master):
        super().__init__(master)

    def create_spectral_plot(self, pc_index=None, d=None):
        """
        Create a spectral plot for a single ROI using the new plotting system.

        Args:
            pc_index: Polygon/ROI index to plot
            d: If not None, plot ratio spectrum (divided by denominator)
        """
        self.num_idx = pc_index
        self.add_spectrum_flag = d is not None

        if pc_index is not None:
            if d is not None:
                self.spectrum = (
                    self.polygon_spectra[pc_index]
                    / self.polygon_spectra[self.denominator_index]
                )
                label_suffix = " (ratio)"
            else:
                self.spectrum = self.polygon_spectra[pc_index]
                label_suffix = ""
            self.spectrum_label = str(self.polygons_table_data[pc_index][0]) + label_suffix

        # Store the original spectrum for correction toggle
        self.original_spectrum = self.spectrum.copy()

        # Use the new plotting system
        self.plot_manager.cleanup_closed_windows()
        active_plot = self.plot_manager.get_active_plot()

        if active_plot is None:
            # Create a new plot window
            active_plot = self._create_spectral_plot_window(title="ROI Spectral Plot")

        # Check accumulate mode
        if not active_plot.is_accumulate_mode():
            active_plot.clear_spectra(update_view=False)

        # Add the spectrum
        spectrum_data = SpectrumData(
            wavelengths=np.array(self.left_wvl),
            values=self.spectrum.copy(),
            label=self.spectrum_label,
            color=self.polygon_colors[pc_index] if pc_index < len(self.polygon_colors) else 'blue',
            metadata={'source': 'roi', 'roi_index': pc_index, 'is_ratio': d is not None}
        )
        active_plot.add_spectrum(spectrum_data)

        # Update backwards-compatibility references
        self.spectral_window = active_plot.window
        self.spectral_ax = active_plot._ax
        self.spectral_canvas = active_plot._canvas


class single_spectrum_window(SpectralCubeAnalysisTool):
    def __init__(self, master):
        super().__init__(master)

    def create_window(self, pc_index):
        self.pc_index = pc_index
        self.single_spectrum_window = tk.Toplevel(self.root)
        self.single_spectrum_window.title("Single Spectrum Plot")

        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.single_spectrum_window)
        ui_frame.pack(side=tk.RIGHT, fill=tk.BOTH)

        spectrum = self.polygon_spectra[self.pc_index]
        single_spectrum_figure, self.single_spectrum_ax = plt.subplots(figsize=(5, 3))
        (self.single_spectrum_line,) = self.single_spectrum_ax.plot(
            self.left_wvl, spectrum, label=self.polygons_table_data[self.pc_index][0]
        )
        self.single_spectrum_ax.legend(loc="best")

        self.single_spectrum_ax.set_xlabel("Wavelength")
        self.single_spectrum_ax.set_ylabel("Value")

        xmin, xmax = self.single_spectrum_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(spectrum[xmin_idx:xmax_idx]), np.nanmax(
            spectrum[xmin_idx:xmax_idx]
        )
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.single_spectrum_ax.set_ylim(min_y - buffer, max_y + buffer)

        self.single_spectral_canvas = FigureCanvasTkAgg(
            single_spectrum_figure, master=self.single_spectrum_window
        )
        self.single_spectral_canvas.get_tk_widget().pack(
            side=tk.TOP, fill=tk.BOTH, expand=True
        )

        # add a button to turn the legend on or off
        self.toggle_legend_button = tk.Button(
            ui_frame, text="Legend (On)", command=self.toggle_spectral_plot_legend
        )
        self.toggle_legend_button.pack(side=tk.TOP)

        # Add a button to reset x-axis span
        self.reset_spectral_x_axis_button = tk.Button(
            ui_frame, text="Reset X-Axis Span", command=self.reset_x_axis_span
        )
        self.reset_spectral_x_axis_button.pack(side=tk.TOP)

        # Add built-in span options to a dropdown menu
        span_options = [
            "Full Span",
            "410 - 1000 nm",
            "410 - 2500 nm",
            "1000 - 2600 nm",
            "1200 - 2000 nm",
            "1800 - 2500 nm",
            "2000 - 2500 nm",
            "2700 - 3900 nm",
        ]
        self.span_var = tk.StringVar()
        self.span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(
            ui_frame, textvariable=self.span_var, values=span_options, state="readonly"
        )
        span_menu.pack(side=tk.TOP)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_x_axis_span)

        # add a button to ignore bad bands
        self.ignore_bad_bands_button = tk.Button(
            ui_frame,
            text="Ignore Bad Bands (Off)",
            command=self.toggle_ignore_bad_bands,
        )
        self.ignore_bad_bands_button.pack(side=tk.TOP)

        # USGS LIBRARY SPECTRA
        # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        self.library_spectra_var = tk.StringVar(ui_frame)
        self.library_spectra_var.set("None")  # default value
        self.library_spectra_list = [
            name.split("/")[-1] for name in self.usgs_spectra_folders
        ]
        # sort library_spectra_list alphabetically
        self.library_spectra_list.sort()
        # label the librar spectra drop down menu "library spectra"
        self.library_spectra_label = tk.Label(ui_frame, text="Library Spectra")
        self.library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.library_spectra_dropdown = tk.OptionMenu(
            ui_frame,
            self.library_spectra_var,
            *self.library_spectra_list,
            command=self.plot_library_spectra,
        )
        self.library_spectra_dropdown.pack(side=tk.TOP)

        # MICA LIBRARY SPECTRA
        # add a drop down menu to plot the library spectra, contents of the drop down menu are the library spectra located in the self.mica_library_spectra folder
        self.mica_library_spectra_var = tk.StringVar(ui_frame)
        self.mica_library_spectra_var.set("None")  # default value
        self.mica_library_spectra_list = [
            name.split(".")[0].split("/")[-1] for name in self.mica_spectra_names
        ]
        # sort mica_library_spectra_list alphabetically
        self.mica_library_spectra_list.sort()
        self.mica_library_spectra_label = tk.Label(
            ui_frame, text="MICA Library Spectra"
        )
        self.mica_library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.mica_library_spectra_dropdown = tk.OptionMenu(
            ui_frame,
            self.mica_library_spectra_var,
            *self.mica_library_spectra_list,
            command=self.plot_mica_library_spectra,
        )
        self.mica_library_spectra_dropdown.pack(side=tk.TOP)

        # add a button to remove library spectra from the plot
        self.remove_library_spectra_button = tk.Button(
            ui_frame, text="Remove Library Spectra", command=self.remove_library_spectra
        )
        self.remove_library_spectra_button.pack(side=tk.TOP)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(
            self.single_spectral_canvas, self.single_spectrum_window
        )
        toolbar.update()
        self.single_spectral_canvas.get_tk_widget().pack(
            side=tk.TOP, fill=tk.BOTH, expand=True
        )

        # Add span selector for x-axis
        self.create_x_axis_span_selector(self.left_wvl)


if __name__ == "__main__":
    root = tk.Tk()
    app = SpectralCubeAnalysisTool(root)
    root.protocol("WM_DELETE_WINDOW", app.on_closing)
    root.mainloop()
