"""
SpectralPlotWindow: Self-contained spectral plot window.

This module provides an independent spectral plot window that manages its own
spectra, supports active/inactive visual states, and can accept dropped spectra.
"""

import tkinter as tk
import warnings
from tkinter import ttk
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import SpanSelector
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Callable
import matplotlib.colors as mcolors


def _safe_nanmedian(values: np.ndarray, default: float = 0.0) -> float:
    """np.nanmedian, but quiet on all-NaN input and returns `default` then.

    Used for stretch/offset pivots where a NaN pivot would poison the
    downstream `pivot + (original - pivot) * stretch + offset` arithmetic.
    """
    arr = np.asarray(values)
    if arr.size == 0 or np.all(np.isnan(arr)):
        return default
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return float(np.nanmedian(arr))


@dataclass
class SpectrumData:
    """
    Data class representing a single spectrum with metadata.

    Attributes:
        wavelengths: Array of wavelength values (nm)
        values: Array of spectral values (reflectance, etc.)
        label: Display label for the spectrum
        color: Line color (matplotlib color string)
        metadata: Additional metadata (coordinates, source, etc.)
        stretch: Vertical stretch factor (1.0 = no change)
        offset: Vertical offset value (0.0 = no change)
    """
    wavelengths: np.ndarray
    values: np.ndarray
    label: str
    color: str = 'blue'
    metadata: Dict[str, Any] = field(default_factory=dict)
    stretch: float = 1.0
    offset: float = 0.0

    def copy(self) -> 'SpectrumData':
        """Create a deep copy of this spectrum data."""
        return SpectrumData(
            wavelengths=self.wavelengths.copy(),
            values=self.values.copy(),
            label=self.label,
            color=self.color,
            metadata=self.metadata.copy(),
            stretch=self.stretch,
            offset=self.offset
        )


class SpectralPlotWindow:
    """
    Self-contained window for displaying spectral plots.

    This class manages its own spectrum storage, plot updates, and UI controls.
    It integrates with PlotManager for active plot tracking and drag-drop support.

    Attributes:
        window: The Tkinter Toplevel window
        plot_id: Unique identifier assigned by PlotManager
        _spectra: List of SpectrumData currently displayed
        _is_active: Whether this is the active plot
    """

    # Class-level color cycle for automatic spectrum coloring
    _color_cycle = list(mcolors.TABLEAU_COLORS.values())
    _next_color_idx = 0

    # Span preset options
    SPAN_OPTIONS = [
        "Full Span",
        "410 - 1000 nm",
        "410 - 2500 nm",
        "1000 - 2600 nm",
        "1200 - 2000 nm",
        "1800 - 2500 nm",
        "2000 - 2500 nm",
        "2700 - 3900 nm",
    ]

    SPAN_RANGES = {
        "410 - 1000 nm": (410, 1000),
        "410 - 2500 nm": (410, 2500),
        "1000 - 2600 nm": (1000, 2600),
        "1200 - 2000 nm": (1200, 2000),
        "1800 - 2500 nm": (1800, 2500),
        "2000 - 2500 nm": (2000, 2500),
        "2700 - 3900 nm": (2700, 3900),
    }

    def __init__(
        self,
        root: tk.Tk,
        title: str = "Spectral Plot",
        on_close_callback: Optional[Callable[['SpectralPlotWindow'], None]] = None,
        on_focus_callback: Optional[Callable[['SpectralPlotWindow'], None]] = None,
        bad_bands: Optional[np.ndarray] = None,
        bad_bands_provider: Optional[Callable[[], Optional[np.ndarray]]] = None,
        library_spectra_folders: Optional[List[str]] = None,
        mica_spectra_names: Optional[List[str]] = None,
        usgs_spectra_path: str = "librarySpectra/",
        mica_spectra_path: str = "mrocr_8001/mrocr_8001/data/",
    ):
        """
        Create a new spectral plot window.

        Args:
            root: Parent Tkinter root window
            title: Window title
            on_close_callback: Called when window is closed
            on_focus_callback: Called when window receives focus
            bad_bands: Array of band quality flags (1=good, 0=bad). Snapshot
                value; use bad_bands_provider if the parent's bad_bands list
                can be reassigned after this window is created.
            bad_bands_provider: Callable returning the current bad_bands.
                When set, takes precedence over the bad_bands snapshot — this
                is the right choice when the parent application rebinds its
                bad_bands attribute (e.g. after Reset Bad Bands or a bad-band
                table edit), so the toggle reflects the live value rather
                than a stale reference captured at construction.
            library_spectra_folders: List of USGS library spectrum folders
            mica_spectra_names: List of MICA library spectrum file names
            usgs_spectra_path: Path to USGS library spectra
            mica_spectra_path: Path to MICA library spectra
        """
        self.root = root
        self.plot_id: Optional[str] = None
        self._spectra: List[SpectrumData] = []
        self._spectrum_lines: Dict[int, Any] = {}  # spectrum index -> Line2D
        self._is_active = False
        self._on_close_callback = on_close_callback
        self._on_focus_callback = on_focus_callback
        self._base_title = title
        self._bad_bands = bad_bands
        self._bad_bands_provider = bad_bands_provider
        self._ignore_bad_bands = False

        # Library spectra paths and names
        self._usgs_spectra_path = usgs_spectra_path
        self._mica_spectra_path = mica_spectra_path
        self._library_spectra_folders = library_spectra_folders or []
        self._mica_spectra_names = mica_spectra_names or []

        # Library spectra tracking: {name: {'line': Line2D, 'original_y': array, ...}}
        self._library_spectra_data: Dict[str, Dict] = {}
        self._library_controls_window: Optional[tk.Toplevel] = None

        # Accumulate mode: when True, new spectra are added without clearing existing ones
        self._accumulate_mode = False

        # Per-spectrum display data: {idx: {'line': Line2D, 'original_y': array, 'pivot': float}}
        self._spectrum_display_data: Dict[int, Dict] = {}
        self._spectrum_controls_window: Optional[tk.Toplevel] = None

        # Drag state tracking
        self._drag_start_pos: Optional[tuple] = None
        self._potential_drag_idx: Optional[int] = None
        self._is_dragging = False
        self._drag_threshold = 10  # pixels needed to start drag

        # Create the window
        self._create_window(title)

    def _create_window(self, title: str) -> None:
        """Create the Toplevel window and all UI components."""
        self.window = tk.Toplevel(self.root)
        self.window.title(title)
        self.window.geometry("700x500")

        # Bind window events
        self.window.protocol("WM_DELETE_WINDOW", self._on_window_close)
        self.window.bind("<FocusIn>", self._on_window_focus)

        # Create main container with border for active indication
        self._border_frame = tk.Frame(
            self.window,
            bd=1,
            relief=tk.SOLID,
            highlightthickness=0
        )
        self._border_frame.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

        # Create UI frame for controls on the right
        self._ui_frame = tk.Frame(self._border_frame)
        self._ui_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=5, pady=5)

        # Create matplotlib figure and canvas
        self._figure = Figure(figsize=(5, 3))
        self._ax = self._figure.add_subplot(111)
        self._ax.set_xlabel('Wavelength (nm)')
        self._ax.set_ylabel('Value')
        self._ax.set_title('Spectral Plot')

        self._canvas = FigureCanvasTkAgg(self._figure, master=self._border_frame)
        self._canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Bind canvas click for activation
        self._canvas.get_tk_widget().bind("<Button-1>", self._on_canvas_click)

        # Connect mouse events for drag-drop (using button_press instead of pick_event)
        self._canvas.mpl_connect('button_press_event', self._on_button_press)
        self._canvas.mpl_connect('button_release_event', self._on_button_release)
        self._canvas.mpl_connect('motion_notify_event', self._on_motion)

        # Create navigation toolbar
        self._toolbar = NavigationToolbar2Tk(self._canvas, self._border_frame)
        self._toolbar.update()

        # Create UI controls
        self._create_controls()

        # Initialize span selector (will be set up when first spectrum is added)
        self._span_selector = None

        # Draw initial empty state
        self._canvas.draw_idle()

    def _create_controls(self) -> None:
        """Create the control buttons and dropdowns."""
        # Accumulate mode toggle button
        self._accumulate_mode_btn = tk.Button(
            self._ui_frame,
            text="Mode: Replace",
            command=self._toggle_accumulate_mode
        )
        self._accumulate_mode_btn.pack(side=tk.TOP, pady=2)

        # Legend toggle button
        self._legend_on = True
        self._toggle_legend_btn = tk.Button(
            self._ui_frame,
            text="Legend (On)",
            command=self._toggle_legend
        )
        self._toggle_legend_btn.pack(side=tk.TOP, pady=2)

        # Reset X-axis span button
        self._reset_span_btn = tk.Button(
            self._ui_frame,
            text="Reset X-Axis Span",
            command=self._reset_x_span
        )
        self._reset_span_btn.pack(side=tk.TOP, pady=2)

        # Span preset dropdown
        self._span_var = tk.StringVar(value="Full Span")
        self._span_combo = ttk.Combobox(
            self._ui_frame,
            textvariable=self._span_var,
            values=self.SPAN_OPTIONS,
            state="readonly",
            width=15
        )
        self._span_combo.pack(side=tk.TOP, pady=2)
        self._span_combo.bind("<<ComboboxSelected>>", self._on_span_selected)

        # Manual X-axis limit entry fields
        xlim_frame = tk.Frame(self._ui_frame)
        xlim_frame.pack(side=tk.TOP, pady=2)
        tk.Label(xlim_frame, text="X Min:").pack(side=tk.LEFT)
        self._xmin_entry = tk.Entry(xlim_frame, width=6)
        self._xmin_entry.pack(side=tk.LEFT, padx=1)
        tk.Label(xlim_frame, text="Max:").pack(side=tk.LEFT)
        self._xmax_entry = tk.Entry(xlim_frame, width=6)
        self._xmax_entry.pack(side=tk.LEFT, padx=1)
        self._apply_xlim_btn = tk.Button(
            self._ui_frame,
            text="Apply X Limits",
            command=self._apply_manual_xlim
        )
        self._apply_xlim_btn.pack(side=tk.TOP, pady=2)

        # Manual Y-axis limit entry fields
        ylim_frame = tk.Frame(self._ui_frame)
        ylim_frame.pack(side=tk.TOP, pady=2)
        tk.Label(ylim_frame, text="Y Min:").pack(side=tk.LEFT)
        self._ymin_entry = tk.Entry(ylim_frame, width=6)
        self._ymin_entry.pack(side=tk.LEFT, padx=1)
        tk.Label(ylim_frame, text="Max:").pack(side=tk.LEFT)
        self._ymax_entry = tk.Entry(ylim_frame, width=6)
        self._ymax_entry.pack(side=tk.LEFT, padx=1)
        self._apply_ylim_btn = tk.Button(
            self._ui_frame,
            text="Apply Y Limits",
            command=self._apply_manual_ylim
        )
        self._apply_ylim_btn.pack(side=tk.TOP, pady=2)

        # Bad bands toggle button
        self._toggle_bad_bands_btn = tk.Button(
            self._ui_frame,
            text="Ignore Bad Bands (Off)",
            command=self._toggle_bad_bands
        )
        self._toggle_bad_bands_btn.pack(side=tk.TOP, pady=2)

        # 1µm offset correction button (for CRISM data)
        self._offset_correction_applied = False
        self._toggle_offset_btn = tk.Button(
            self._ui_frame,
            text="Correct 1µm Offset (Off)",
            command=self._toggle_offset_correction
        )
        self._toggle_offset_btn.pack(side=tk.TOP, pady=2)

        # Separator
        ttk.Separator(self._ui_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=5)

        # USGS Library Spectra
        tk.Label(self._ui_frame, text="USGS Library Spectra").pack(side=tk.TOP, pady=2)
        self._library_spectra_var = tk.StringVar(value="None")
        library_list = sorted([name.split('/')[-1] for name in self._library_spectra_folders])
        if library_list:
            self._library_dropdown = tk.OptionMenu(
                self._ui_frame,
                self._library_spectra_var,
                *library_list,
                command=self._plot_usgs_library_spectrum
            )
            self._library_dropdown.pack(side=tk.TOP, pady=2)

        # MICA Library Spectra
        tk.Label(self._ui_frame, text="MICA Library Spectra").pack(side=tk.TOP, pady=2)
        self._mica_library_var = tk.StringVar(value="None")
        mica_list = sorted([name.split(".")[0].split("/")[-1] for name in self._mica_spectra_names])
        if mica_list:
            self._mica_dropdown = tk.OptionMenu(
                self._ui_frame,
                self._mica_library_var,
                *mica_list,
                command=self._plot_mica_library_spectrum
            )
            self._mica_dropdown.pack(side=tk.TOP, pady=2)

        # Remove library spectra button
        self._remove_library_btn = tk.Button(
            self._ui_frame,
            text="Remove Library Spectra",
            command=self._remove_library_spectra
        )
        self._remove_library_btn.pack(side=tk.TOP, pady=2)

        # Separator
        ttk.Separator(self._ui_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=5)

        # Spectrum Controls button
        self._spectrum_controls_btn = tk.Button(
            self._ui_frame,
            text="Spectrum Controls",
            command=self._show_spectrum_controls
        )
        self._spectrum_controls_btn.pack(side=tk.TOP, pady=2)

        # Clear all spectra button
        self._clear_btn = tk.Button(
            self._ui_frame,
            text="Clear All Spectra",
            command=self.clear_spectra
        )
        self._clear_btn.pack(side=tk.TOP, pady=2)

    def set_plot_id(self, plot_id: str) -> None:
        """Set the plot ID (called by PlotManager)."""
        self.plot_id = plot_id

    def add_spectrum(self, spectrum: SpectrumData, update_view: bool = True) -> int:
        """
        Add a spectrum to this plot.

        Args:
            spectrum: The SpectrumData to add
            update_view: Whether to update the plot immediately

        Returns:
            Index of the added spectrum
        """
        # Assign color if not specified or if 'blue' (default)
        if spectrum.color == 'blue':
            spectrum.color = self._get_next_color()

        # Apply corrections to displayed values
        values = spectrum.values.copy()

        # Apply offset correction if enabled
        if self._offset_correction_applied:
            values = self._apply_1um_offset_correction(values, spectrum.wavelengths)

        # Apply bad bands masking if enabled
        bb = self._get_bad_bands()
        if self._ignore_bad_bands and bb is not None:
            values = np.where(bb == 0, np.nan, values)

        # Add to internal list
        idx = len(self._spectra)
        self._spectra.append(spectrum)

        # Calculate pivot point for stretch/offset operations.
        # An all-NaN spectrum is legitimate (e.g. a column-locked denominator
        # polygon that landed in a no-data area of a CRISM tile). nanmedian
        # would emit "All-NaN slice encountered" and return NaN, which would
        # then poison stretch/offset arithmetic; default to 0.0 instead.
        pivot = _safe_nanmedian(values)

        # Plot the spectrum
        line, = self._ax.plot(
            spectrum.wavelengths,
            values,
            color=spectrum.color,
            label=spectrum.label,
            picker=5  # Enable picking for drag-drop
        )
        self._spectrum_lines[idx] = line

        # Store display data for stretch/offset controls
        self._spectrum_display_data[idx] = {
            'line': line,
            'original_y': values.copy(),
            'stretch': spectrum.stretch,
            'offset': spectrum.offset,
            'pivot': pivot
        }

        # Set up span selector if this is the first spectrum
        if idx == 0:
            self._setup_span_selector(spectrum.wavelengths)

        if update_view:
            self._update_view()

        return idx

    def remove_spectrum(self, index: int) -> bool:
        """
        Remove a spectrum by index.

        Args:
            index: Index of the spectrum to remove

        Returns:
            True if removed, False if index invalid
        """
        if index < 0 or index >= len(self._spectra):
            return False

        # Remove the line from the plot
        if index in self._spectrum_lines:
            self._spectrum_lines[index].remove()
            del self._spectrum_lines[index]

        # Remove from display data
        if index in self._spectrum_display_data:
            del self._spectrum_display_data[index]

        # Remove from spectra list
        del self._spectra[index]

        # Re-index remaining lines and display data
        new_lines = {}
        new_display_data = {}
        for old_idx, line in self._spectrum_lines.items():
            new_idx = old_idx if old_idx < index else old_idx - 1
            new_lines[new_idx] = line
        for old_idx, data in self._spectrum_display_data.items():
            new_idx = old_idx if old_idx < index else old_idx - 1
            new_display_data[new_idx] = data
        self._spectrum_lines = new_lines
        self._spectrum_display_data = new_display_data

        self._update_view()
        return True

    def clear_spectra(self, update_view: bool = True) -> None:
        """Remove all spectra from the plot.

        Args:
            update_view: Whether to update the plot immediately. Set to False
                        if you plan to add spectra right after clearing.
        """
        for line in self._spectrum_lines.values():
            line.remove()
        self._spectrum_lines.clear()
        self._spectrum_display_data.clear()
        self._spectra.clear()
        if update_view:
            self._update_view()

    def get_spectrum(self, index: int) -> Optional[SpectrumData]:
        """Get a spectrum by index."""
        if 0 <= index < len(self._spectra):
            return self._spectra[index]
        return None

    def get_all_spectra(self) -> List[SpectrumData]:
        """Get all spectra in this plot."""
        return self._spectra.copy()

    def get_spectrum_count(self) -> int:
        """Get the number of spectra in this plot."""
        return len(self._spectra)

    def set_active(self, active: bool) -> None:
        """
        Set the active state of this plot window.

        Args:
            active: True to mark as active, False for inactive
        """
        self._is_active = active
        self._update_active_appearance()

    def is_active(self) -> bool:
        """Check if this plot is the active plot."""
        return self._is_active

    def is_window_valid(self) -> bool:
        """Check if the window still exists."""
        try:
            return self.window.winfo_exists()
        except tk.TclError:
            return False

    def focus(self) -> None:
        """Bring the window to focus."""
        if self.is_window_valid():
            self.window.focus_force()
            self.window.lift()

    def accept_drop(self, spectrum: SpectrumData) -> bool:
        """
        Accept a dropped spectrum from drag-and-drop.

        Args:
            spectrum: The spectrum being dropped

        Returns:
            True if accepted, False otherwise
        """
        # Make a copy to avoid sharing state
        self.add_spectrum(spectrum.copy())
        return True

    def _update_active_appearance(self) -> None:
        """Update the visual appearance based on active state."""
        if self._is_active:
            self._border_frame.config(
                highlightbackground='cyan',
                highlightcolor='cyan',
                highlightthickness=3
            )
            self.window.title(f"[ACTIVE] {self._base_title}")
        else:
            self._border_frame.config(
                highlightbackground='gray',
                highlightcolor='gray',
                highlightthickness=1
            )
            self.window.title(self._base_title)

    def _update_view(self, immediate: bool = False) -> None:
        """Update the plot view (legend, limits, redraw).

        Args:
            immediate: If True, use draw() for immediate update.
                      If False (default), use draw_idle() for better performance.
        """
        # Update legend - use 'upper right' for faster rendering than 'best'
        if self._legend_on and self._spectra:
            self._ax.legend(loc='upper right')
        elif self._ax.get_legend():
            self._ax.get_legend().remove()

        # Auto-scale y-axis based on visible data
        self._auto_scale_y()

        # Use draw_idle for non-blocking updates (faster)
        if immediate:
            self._canvas.draw_idle()
        else:
            self._canvas.draw_idle()

    def _auto_scale_y(self) -> None:
        """Auto-scale the y-axis based on visible x range."""
        if not self._spectra:
            return

        xmin, xmax = self._ax.get_xlim()
        all_min_y = []
        all_max_y = []

        for spectrum in self._spectra:
            # Find indices within x range
            wvl = spectrum.wavelengths
            mask = (wvl >= xmin) & (wvl <= xmax)
            if not np.any(mask):
                continue

            values = spectrum.values[mask]
            bb = self._get_bad_bands()
            if self._ignore_bad_bands and bb is not None:
                # Get corresponding bad band mask
                bb_mask = bb[mask] if len(bb) == len(spectrum.values) else None
                if bb_mask is not None:
                    values = np.where(bb_mask == 0, np.nan, values)

            if not np.all(np.isnan(values)):
                all_min_y.append(np.nanmin(values))
                all_max_y.append(np.nanmax(values))

        if all_min_y and all_max_y:
            min_y = min(all_min_y)
            max_y = max(all_max_y)
            buffer = (max_y - min_y) * 0.1
            self._ax.set_ylim(min_y - buffer, max_y + buffer)

    def _setup_span_selector(self, wavelengths: np.ndarray) -> None:
        """Set up the span selector for x-axis zoom."""
        self._wavelengths = wavelengths

        def on_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(wavelengths - xmin))
            xmax_idx = np.argmin(np.abs(wavelengths - xmax))
            self._ax.set_xlim(wavelengths[xmin_idx], wavelengths[xmax_idx])
            self._auto_scale_y()
            self._canvas.draw_idle()

        self._span_selector = SpanSelector(
            self._ax,
            on_span_select,
            "horizontal",
            useblit=True
        )

    def _get_next_color(self) -> str:
        """Get the next color from the color cycle."""
        color = self._color_cycle[SpectralPlotWindow._next_color_idx % len(self._color_cycle)]
        SpectralPlotWindow._next_color_idx += 1
        return color

    # ----- Event handlers -----

    def _on_window_close(self) -> None:
        """Handle window close event."""
        # Close library controls window if open
        if self._library_controls_window and self._library_controls_window.winfo_exists():
            self._library_controls_window.destroy()

        # Close spectrum controls window if open
        if self._spectrum_controls_window and self._spectrum_controls_window.winfo_exists():
            self._spectrum_controls_window.destroy()

        if self._on_close_callback:
            self._on_close_callback(self)
        self.window.destroy()

    def _on_window_focus(self, event=None) -> None:
        """Handle window focus event."""
        if self._on_focus_callback:
            self._on_focus_callback(self)

    def _on_canvas_click(self, event=None) -> None:
        """Handle canvas click for activation."""
        if self._on_focus_callback:
            self._on_focus_callback(self)

    def _find_spectrum_at_point(self, event) -> Optional[int]:
        """Find which spectrum (if any) is under the mouse cursor."""
        if event.inaxes != self._ax:
            return None

        for idx, line in self._spectrum_lines.items():
            contains, _ = line.contains(event)
            if contains:
                return idx
        return None

    def _on_button_press(self, event) -> None:
        """Handle mouse button press for potential drag start."""
        if event.button != 1:  # Only left button
            return

        idx = self._find_spectrum_at_point(event)
        if idx is not None:
            self._potential_drag_idx = idx
            # Store pixel coordinates for threshold calculation
            self._drag_start_pos = (event.x, event.y)
        else:
            self._potential_drag_idx = None
            self._drag_start_pos = None

    def _on_motion(self, event) -> None:
        """Handle mouse motion for drag detection."""
        if self._drag_start_pos is None or self._is_dragging:
            return

        if event.x is None or event.y is None:
            return

        # Calculate distance moved
        dx = abs(event.x - self._drag_start_pos[0])
        dy = abs(event.y - self._drag_start_pos[1])

        # Check if we've moved past the threshold
        if dx > self._drag_threshold or dy > self._drag_threshold:
            self._is_dragging = True
            self._start_drag_operation(event)

    def _on_button_release(self, event) -> None:
        """Handle mouse button release."""
        if event.button != 1:
            return

        # Reset drag state
        self._drag_start_pos = None
        self._potential_drag_idx = None
        self._is_dragging = False

    def _start_drag_operation(self, event) -> None:
        """Start the actual drag operation once threshold is exceeded."""
        if self._potential_drag_idx is None or self._potential_drag_idx >= len(self._spectra):
            return

        spectrum = self._spectra[self._potential_drag_idx]

        # Get screen coordinates
        canvas_widget = self._canvas.get_tk_widget()
        screen_x = canvas_widget.winfo_rootx() + int(event.x)
        screen_y = canvas_widget.winfo_rooty() + int(event.y)

        try:
            from drag_drop_manager import get_drag_drop_manager
            drag_manager = get_drag_drop_manager()
            drag_manager.start_drag(self, spectrum, screen_x, screen_y)
        except Exception as e:
            print(f"Error starting drag: {e}")

    # ----- Control handlers -----

    def _toggle_accumulate_mode(self) -> None:
        """Toggle between replace and accumulate modes."""
        self._accumulate_mode = not self._accumulate_mode
        if self._accumulate_mode:
            self._accumulate_mode_btn.config(text="Mode: Accumulate", relief=tk.SUNKEN)
        else:
            self._accumulate_mode_btn.config(text="Mode: Replace", relief=tk.RAISED)

    def is_accumulate_mode(self) -> bool:
        """Check if accumulate mode is enabled."""
        return self._accumulate_mode

    def _toggle_legend(self) -> None:
        """Toggle the legend visibility."""
        self._legend_on = not self._legend_on
        if self._legend_on:
            self._toggle_legend_btn.config(text="Legend (On)", relief=tk.RAISED)
        else:
            self._toggle_legend_btn.config(text="Legend (Off)", relief=tk.SUNKEN)
        self._update_view()

    def _reset_x_span(self) -> None:
        """Reset x-axis to full span."""
        if self._spectra:
            wvl = self._spectra[0].wavelengths
            self._ax.set_xlim(wvl[0], wvl[-1])
            self._auto_scale_y()
            self._canvas.draw_idle()
            self._span_var.set("Full Span")

    def _on_span_selected(self, event=None) -> None:
        """Handle span preset selection."""
        selected = self._span_var.get()
        if selected == "Full Span":
            self._reset_x_span()
        elif selected in self.SPAN_RANGES:
            xmin, xmax = self.SPAN_RANGES[selected]
            self._ax.set_xlim(xmin, xmax)
            self._auto_scale_y()
            self._canvas.draw_idle()

    def _apply_manual_xlim(self) -> None:
        """Apply manually entered X-axis limits."""
        try:
            xmin_str = self._xmin_entry.get().strip()
            xmax_str = self._xmax_entry.get().strip()
            if not xmin_str or not xmax_str:
                return
            xmin = float(xmin_str)
            xmax = float(xmax_str)
            if xmin >= xmax:
                return
            self._ax.set_xlim(xmin, xmax)
            self._auto_scale_y()
            self._canvas.draw_idle()
            self._span_var.set("Full Span")  # Reset dropdown since custom range
        except ValueError:
            pass  # Invalid input, ignore

    def _apply_manual_ylim(self) -> None:
        """Apply manually entered Y-axis limits."""
        try:
            ymin_str = self._ymin_entry.get().strip()
            ymax_str = self._ymax_entry.get().strip()
            if not ymin_str or not ymax_str:
                return
            ymin = float(ymin_str)
            ymax = float(ymax_str)
            if ymin >= ymax:
                return
            self._ax.set_ylim(ymin, ymax)
            self._canvas.draw_idle()
        except ValueError:
            pass  # Invalid input, ignore

    def _get_bad_bands(self) -> Optional[np.ndarray]:
        """Resolve the current bad-bands flags as an ndarray, or None.

        Prefers the provider callback (which returns the parent app's
        live attribute) over the construction-time snapshot, so this
        survives the parent rebinding its bad_bands list (Reset Bad
        Bands, table edits, etc.).
        """
        bb = self._bad_bands_provider() if self._bad_bands_provider else self._bad_bands
        if bb is None:
            return None
        return np.asarray(bb)

    def _toggle_bad_bands(self) -> None:
        """Toggle bad bands masking."""
        self._ignore_bad_bands = not self._ignore_bad_bands

        if self._ignore_bad_bands:
            self._toggle_bad_bands_btn.config(text="Ignore Bad Bands (On)", relief=tk.SUNKEN)
        else:
            self._toggle_bad_bands_btn.config(text="Ignore Bad Bands (Off)", relief=tk.RAISED)

        # Replot all spectra with new masking
        self._replot_all_spectra()

    def _toggle_offset_correction(self) -> None:
        """Toggle the 1µm offset correction for CRISM spectra."""
        self._offset_correction_applied = not self._offset_correction_applied

        if self._offset_correction_applied:
            self._toggle_offset_btn.config(text="Correct 1µm Offset (On)", relief=tk.SUNKEN)
        else:
            self._toggle_offset_btn.config(text="Correct 1µm Offset (Off)", relief=tk.RAISED)

        # Replot all spectra with correction
        self._replot_all_spectra()

    def _apply_1um_offset_correction(self, spectrum: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
        """
        Correct the 1µm offset in CRISM spectra.

        This corrects the discontinuity between short-wave and long-wave detectors
        around 1000nm by adjusting the short-wavelength values to align with the
        long-wavelength values.
        """
        # Find the wavelength closest to 1000nm (detector boundary)
        target_wavelength = 1000
        detector_boundary_idx = np.argmin(np.abs(wavelengths - target_wavelength))

        # Define windows around the detector boundary (5 samples each side)
        window_size = 5

        sw_start = max(0, detector_boundary_idx - window_size)
        sw_end = detector_boundary_idx
        lw_start = detector_boundary_idx + 1
        lw_end = min(len(spectrum), detector_boundary_idx + window_size + 1)

        sw_values = spectrum[sw_start:sw_end]
        lw_values = spectrum[lw_start:lw_end]

        if len(sw_values) > 0 and len(lw_values) > 0:
            sw_median = np.nanmedian(sw_values)
            lw_median = np.nanmedian(lw_values)

            if not np.isnan(sw_median) and not np.isnan(lw_median):
                offset = lw_median - sw_median
                corrected = spectrum.copy()
                corrected[:detector_boundary_idx] += offset
                return corrected

        return spectrum

    def _replot_all_spectra(self) -> None:
        """Replot all spectra (used after bad bands or offset correction toggle)."""
        bb = self._get_bad_bands()
        for idx, spectrum in enumerate(self._spectra):
            values = spectrum.values.copy()

            # Apply offset correction if enabled
            if self._offset_correction_applied:
                values = self._apply_1um_offset_correction(values, spectrum.wavelengths)

            # Apply bad bands masking if enabled
            if self._ignore_bad_bands and bb is not None:
                values = np.where(bb == 0, np.nan, values)

            if idx in self._spectrum_lines:
                self._spectrum_lines[idx].set_ydata(values)

        self._update_view()

    def _plot_usgs_library_spectrum(self, spec_name: str) -> None:
        """Plot a USGS library spectrum."""
        try:
            from hypyrameter.utils import getWavelengthFromUSGS, getReflectanceFromUSGS

            path_to_spec = f"{self._usgs_spectra_path}{spec_name}/{spec_name}.txt"
            path_to_wvl = f"{self._usgs_spectra_path}{spec_name}/*Wavelengths*.txt"

            library_wvl = getWavelengthFromUSGS(path_to_wvl)
            library_refl = getReflectanceFromUSGS(path_to_spec)

            self._add_library_spectrum(spec_name, library_wvl, library_refl)
        except Exception as e:
            print(f"Error loading USGS library spectrum: {e}")

    def _plot_mica_library_spectrum(self, spec_name: str) -> None:
        """Plot a MICA library spectrum."""
        try:
            import pandas as pd

            spec_path = f"{self._mica_spectra_path}{spec_name}.tab"
            # MICA files use comma separators, not tabs
            data = pd.read_csv(spec_path, sep=',', header=None)

            # MICA format: column 0 is wavelength (microns), column 1 is reflectance
            library_wvl = data[0].values * 1000  # Convert to nm
            library_refl = data[1].values

            self._add_library_spectrum(spec_name, library_wvl, library_refl)
        except Exception as e:
            print(f"Error loading MICA library spectrum: {e}")

    def _add_library_spectrum(self, name: str, wvl: np.ndarray, refl: np.ndarray) -> None:
        """Add a library spectrum to the plot with scaling."""
        if not self._spectra:
            print("No data spectra to scale library spectrum against")
            return

        # Get current view limits for scaling
        xmin, xmax = self._ax.get_xlim()

        # Get the first spectrum's values in the visible range for scaling
        ref_spectrum = self._spectra[0]
        ref_wvl = ref_spectrum.wavelengths
        xmin_idx = np.argmin(np.abs(ref_wvl - xmin))
        xmax_idx = np.argmin(np.abs(ref_wvl - xmax))
        ref_slice = ref_spectrum.values[xmin_idx:xmax_idx]

        if np.all(np.isnan(ref_slice)):
            print("Cannot scale library spectrum: reference data is all NaN")
            return

        # Scale library spectrum to match data range
        scaled_refl = (refl - np.nanmin(refl)) * (
            (np.nanmax(ref_slice) - np.nanmin(ref_slice)) /
            (np.nanmax(refl) - np.nanmin(refl))
        ) + np.nanmin(ref_slice)

        # Plot the library spectrum
        line, = self._ax.plot(wvl, scaled_refl, label=name, linestyle='--')

        # Store for controls
        pivot = np.nanmedian(scaled_refl)
        self._library_spectra_data[name] = {
            'line': line,
            'original_y': scaled_refl.copy(),
            'wvl': wvl,
            'stretch': 1.0,
            'offset': 0.0,
            'anchor': pivot,
            'pivot': pivot
        }

        self._update_view()
        self._update_library_controls()

    def _remove_library_spectra(self) -> None:
        """Remove all library spectra from the plot."""
        for name, data in self._library_spectra_data.items():
            data['line'].remove()
        self._library_spectra_data.clear()

        if self._library_controls_window and self._library_controls_window.winfo_exists():
            self._library_controls_window.destroy()
            self._library_controls_window = None

        self._update_view()

    def _update_library_controls(self) -> None:
        """Update the library spectra stretch/offset controls window."""
        if not self._library_spectra_data:
            if self._library_controls_window and self._library_controls_window.winfo_exists():
                self._library_controls_window.destroy()
                self._library_controls_window = None
            return

        # Create controls window if needed
        if self._library_controls_window is None or not self._library_controls_window.winfo_exists():
            self._library_controls_window = tk.Toplevel(self.window)
            self._library_controls_window.title("Library Spectra Controls")
            self._library_controls_window.geometry("400x300")

        # Clear existing controls
        for widget in self._library_controls_window.winfo_children():
            widget.destroy()

        # Create scrollable frame
        canvas = tk.Canvas(self._library_controls_window)
        scrollbar = ttk.Scrollbar(self._library_controls_window, orient="vertical", command=canvas.yview)
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
        for spec_name, spec_data in self._library_spectra_data.items():
            frame = ttk.LabelFrame(scrollable_frame, text=spec_name, padding=5)
            frame.pack(fill=tk.X, padx=5, pady=5)

            # Stretch control
            stretch_frame = ttk.Frame(frame)
            stretch_frame.pack(fill=tk.X)
            ttk.Label(stretch_frame, text="Stretch:").pack(side=tk.LEFT)
            stretch_var = tk.DoubleVar(value=spec_data['stretch'])
            stretch_scale = ttk.Scale(
                stretch_frame,
                from_=0.1, to=3.0,
                variable=stretch_var,
                orient=tk.HORIZONTAL,
                command=lambda v, n=spec_name, sv=stretch_var: self._on_library_stretch_change(n, sv.get())
            )
            stretch_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

            # Offset control
            offset_frame = ttk.Frame(frame)
            offset_frame.pack(fill=tk.X)
            ttk.Label(offset_frame, text="Offset:").pack(side=tk.LEFT)
            offset_var = tk.DoubleVar(value=spec_data['offset'])
            offset_scale = ttk.Scale(
                offset_frame,
                from_=-0.5, to=0.5,
                variable=offset_var,
                orient=tk.HORIZONTAL,
                command=lambda v, n=spec_name, ov=offset_var: self._on_library_offset_change(n, ov.get())
            )
            offset_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

    def _on_library_stretch_change(self, name: str, stretch: float) -> None:
        """Handle library spectrum stretch change."""
        if name not in self._library_spectra_data:
            return

        spec_data = self._library_spectra_data[name]
        spec_data['stretch'] = stretch

        # Apply stretch around pivot point
        pivot = spec_data['pivot']
        original = spec_data['original_y']
        new_y = pivot + (original - pivot) * stretch + spec_data['offset']

        spec_data['line'].set_ydata(new_y)
        self._canvas.draw_idle()

    def _on_library_offset_change(self, name: str, offset: float) -> None:
        """Handle library spectrum offset change."""
        if name not in self._library_spectra_data:
            return

        spec_data = self._library_spectra_data[name]
        spec_data['offset'] = offset

        # Apply offset with current stretch
        pivot = spec_data['pivot']
        original = spec_data['original_y']
        new_y = pivot + (original - pivot) * spec_data['stretch'] + offset

        spec_data['line'].set_ydata(new_y)
        self._canvas.draw_idle()

    # ----- Spectrum Controls -----

    def _show_spectrum_controls(self) -> None:
        """Show or update the spectrum controls window."""
        if not self._spectra:
            return

        # Create controls window if needed
        if self._spectrum_controls_window is None or not self._spectrum_controls_window.winfo_exists():
            self._spectrum_controls_window = tk.Toplevel(self.window)
            self._spectrum_controls_window.title("Spectrum Controls")
            self._spectrum_controls_window.geometry("400x400")

        self._update_spectrum_controls()

    def _update_spectrum_controls(self) -> None:
        """Update the spectrum controls window with current spectra."""
        if self._spectrum_controls_window is None or not self._spectrum_controls_window.winfo_exists():
            return

        # Clear existing controls
        for widget in self._spectrum_controls_window.winfo_children():
            widget.destroy()

        if not self._spectra:
            tk.Label(self._spectrum_controls_window, text="No spectra loaded").pack(pady=20)
            return

        # Create scrollable frame
        canvas = tk.Canvas(self._spectrum_controls_window)
        scrollbar = ttk.Scrollbar(self._spectrum_controls_window, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Create controls for each spectrum
        for idx, spectrum in enumerate(self._spectra):
            if idx not in self._spectrum_display_data:
                continue

            spec_data = self._spectrum_display_data[idx]

            frame = ttk.LabelFrame(scrollable_frame, text=spectrum.label, padding=5)
            frame.pack(fill=tk.X, padx=5, pady=5)

            # Color indicator and picker button
            color_frame = ttk.Frame(frame)
            color_frame.pack(fill=tk.X)
            ttk.Label(color_frame, text="Color:").pack(side=tk.LEFT)
            color_btn = tk.Button(
                color_frame,
                text="    ",
                bg=spectrum.color,
                command=lambda i=idx: self._pick_spectrum_color(i)
            )
            color_btn.pack(side=tk.LEFT, padx=5)

            # Stretch control
            stretch_frame = ttk.Frame(frame)
            stretch_frame.pack(fill=tk.X)
            ttk.Label(stretch_frame, text="Stretch:").pack(side=tk.LEFT)
            stretch_var = tk.DoubleVar(value=spec_data['stretch'])
            stretch_scale = ttk.Scale(
                stretch_frame,
                from_=0.1, to=3.0,
                variable=stretch_var,
                orient=tk.HORIZONTAL,
                command=lambda v, i=idx, sv=stretch_var: self._on_spectrum_stretch_change(i, sv.get())
            )
            stretch_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

            # Offset control
            offset_frame = ttk.Frame(frame)
            offset_frame.pack(fill=tk.X)
            ttk.Label(offset_frame, text="Offset:").pack(side=tk.LEFT)
            offset_var = tk.DoubleVar(value=spec_data['offset'])
            offset_scale = ttk.Scale(
                offset_frame,
                from_=-0.5, to=0.5,
                variable=offset_var,
                orient=tk.HORIZONTAL,
                command=lambda v, i=idx, ov=offset_var: self._on_spectrum_offset_change(i, ov.get())
            )
            offset_scale.pack(side=tk.LEFT, fill=tk.X, expand=True)

            # Remove button
            remove_btn = tk.Button(
                frame,
                text="Remove",
                command=lambda i=idx: self._remove_spectrum_and_update(i)
            )
            remove_btn.pack(side=tk.RIGHT, pady=2)

    def _on_spectrum_stretch_change(self, idx: int, stretch: float) -> None:
        """Handle spectrum stretch change."""
        if idx not in self._spectrum_display_data:
            return

        spec_data = self._spectrum_display_data[idx]
        spec_data['stretch'] = stretch

        # Update the SpectrumData object as well
        if idx < len(self._spectra):
            self._spectra[idx].stretch = stretch

        # Apply stretch around pivot point
        pivot = spec_data['pivot']
        original = spec_data['original_y']
        new_y = pivot + (original - pivot) * stretch + spec_data['offset']

        spec_data['line'].set_ydata(new_y)
        self._canvas.draw_idle()

    def _on_spectrum_offset_change(self, idx: int, offset: float) -> None:
        """Handle spectrum offset change."""
        if idx not in self._spectrum_display_data:
            return

        spec_data = self._spectrum_display_data[idx]
        spec_data['offset'] = offset

        # Update the SpectrumData object as well
        if idx < len(self._spectra):
            self._spectra[idx].offset = offset

        # Apply offset with current stretch
        pivot = spec_data['pivot']
        original = spec_data['original_y']
        new_y = pivot + (original - pivot) * spec_data['stretch'] + offset

        spec_data['line'].set_ydata(new_y)
        self._canvas.draw_idle()

    def _pick_spectrum_color(self, idx: int) -> None:
        """Open color picker for a spectrum."""
        from tkinter import colorchooser

        if idx >= len(self._spectra):
            return

        current_color = self._spectra[idx].color
        color = colorchooser.askcolor(initialcolor=current_color, title="Choose Spectrum Color")

        if color[1]:  # color[1] is the hex string
            self._spectra[idx].color = color[1]
            if idx in self._spectrum_display_data:
                self._spectrum_display_data[idx]['line'].set_color(color[1])
            self._canvas.draw_idle()
            # Refresh controls window to update color button
            self._update_spectrum_controls()

    def _remove_spectrum_and_update(self, idx: int) -> None:
        """Remove a spectrum and update the controls window."""
        self.remove_spectrum(idx)
        self._update_spectrum_controls()
