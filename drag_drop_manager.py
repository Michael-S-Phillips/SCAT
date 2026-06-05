"""
DragDropManager: Coordinates drag-and-drop operations between plot windows.

This module handles the visual feedback and coordination for dragging spectra
from one plot window to another.
"""

import tkinter as tk
from typing import Optional, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from spectral_plot_window import SpectralPlotWindow, SpectrumData


class DragDropManager:
    """
    Manages drag-and-drop operations for spectra between plot windows.

    Responsibilities:
    - Track drag state (source plot, spectrum, position)
    - Display floating indicator during drag
    - Highlight valid drop targets
    - Execute drop operation (copy spectrum to target)
    """

    _instance = None

    def __new__(cls, root: tk.Tk = None):
        """Singleton pattern."""
        if cls._instance is None:
            if root is None:
                raise ValueError("Root window required for first instantiation")
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self, root: tk.Tk = None):
        if self._initialized:
            return
        self._initialized = True

        self._root = root
        self._is_dragging = False
        self._source_plot: Optional['SpectralPlotWindow'] = None
        self._dragged_spectrum: Optional['SpectrumData'] = None
        self._drag_indicator: Optional[tk.Toplevel] = None
        self._start_x = 0
        self._start_y = 0

        # Import here to avoid circular imports
        from plot_manager import get_plot_manager
        self._plot_manager = get_plot_manager

    def start_drag(
        self,
        source_plot: 'SpectralPlotWindow',
        spectrum: 'SpectrumData',
        screen_x: int,
        screen_y: int
    ) -> None:
        """
        Begin a drag operation.

        Args:
            source_plot: The plot window where the drag started
            spectrum: The spectrum being dragged
            screen_x: Screen X coordinate of drag start
            screen_y: Screen Y coordinate of drag start
        """
        if self._is_dragging:
            self.cancel_drag()

        self._is_dragging = True
        self._source_plot = source_plot
        self._dragged_spectrum = spectrum
        self._start_x = screen_x
        self._start_y = screen_y

        # Create the drag indicator
        self._create_drag_indicator(spectrum, screen_x, screen_y)

        # Bind mouse events to root and all plot windows
        # This ensures we capture events even when dragging over different windows
        self._root.bind("<B1-Motion>", self._on_drag_motion)
        self._root.bind("<ButtonRelease-1>", self._on_drag_release)

        # Also bind to all plot windows to capture events on Toplevel windows
        plot_manager = self._plot_manager()
        for plot in plot_manager.get_all_plots().values():
            if plot.is_window_valid():
                try:
                    plot.window.bind("<B1-Motion>", self._on_drag_motion)
                    plot.window.bind("<ButtonRelease-1>", self._on_drag_release)
                except tk.TclError:
                    pass

    def _create_drag_indicator(
        self,
        spectrum: 'SpectrumData',
        screen_x: int,
        screen_y: int
    ) -> None:
        """Create a floating indicator showing the dragged spectrum."""
        self._drag_indicator = tk.Toplevel(self._root)
        self._drag_indicator.overrideredirect(True)  # No window decorations
        self._drag_indicator.attributes('-topmost', True)
        self._drag_indicator.attributes('-alpha', 0.8)

        # Create a small preview
        frame = tk.Frame(
            self._drag_indicator,
            bg='white',
            bd=2,
            relief=tk.RAISED
        )
        frame.pack(fill=tk.BOTH, expand=True)

        # Mini matplotlib preview
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

            fig = Figure(figsize=(1.5, 1), dpi=50)
            ax = fig.add_subplot(111)
            ax.plot(
                spectrum.wavelengths,
                spectrum.values,
                color=spectrum.color,
                linewidth=1
            )
            ax.set_xticks([])
            ax.set_yticks([])
            fig.tight_layout(pad=0.1)

            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.draw()
            canvas.get_tk_widget().pack()

        except Exception:
            # Fallback to simple label if matplotlib fails
            label = tk.Label(
                frame,
                text=f"📊 {spectrum.label}",
                bg='white',
                padx=5,
                pady=5
            )
            label.pack()

        # Position the indicator
        self._drag_indicator.geometry(f"+{screen_x + 10}+{screen_y + 10}")

    def _on_drag_motion(self, event) -> None:
        """Handle mouse motion during drag."""
        if not self._is_dragging or not self._drag_indicator:
            return

        # Move the indicator
        screen_x = event.x_root
        screen_y = event.y_root
        self._drag_indicator.geometry(f"+{screen_x + 10}+{screen_y + 10}")

        # Highlight drop targets
        self._update_drop_target_highlights(screen_x, screen_y)

    def _on_drag_release(self, event) -> None:
        """Handle mouse release to complete or cancel drop."""
        if not self._is_dragging:
            return

        screen_x = event.x_root
        screen_y = event.y_root

        # Find the target plot under the cursor
        target_plot = self._find_plot_at_position(screen_x, screen_y)

        if target_plot and target_plot != self._source_plot:
            # Execute the drop
            self._execute_drop(target_plot)

        # Clean up
        self.cancel_drag()

    def _find_plot_at_position(self, screen_x: int, screen_y: int) -> Optional['SpectralPlotWindow']:
        """Find the plot window at the given screen position."""
        plot_manager = self._plot_manager()

        for plot_id, plot in plot_manager.get_all_plots().items():
            if not plot.is_window_valid():
                continue

            try:
                win = plot.window
                x = win.winfo_rootx()
                y = win.winfo_rooty()
                w = win.winfo_width()
                h = win.winfo_height()

                if x <= screen_x <= x + w and y <= screen_y <= y + h:
                    return plot
            except tk.TclError:
                continue

        return None

    def _update_drop_target_highlights(self, screen_x: int, screen_y: int) -> None:
        """Update visual highlights on potential drop targets."""
        plot_manager = self._plot_manager()
        target = self._find_plot_at_position(screen_x, screen_y)

        for plot_id, plot in plot_manager.get_all_plots().items():
            if not plot.is_window_valid():
                continue

            if plot == target and plot != self._source_plot:
                # Highlight as valid drop target
                try:
                    plot._border_frame.config(
                        highlightbackground='green',
                        highlightcolor='green',
                        highlightthickness=4
                    )
                except (tk.TclError, AttributeError):
                    pass
            else:
                # Restore normal appearance
                plot._update_active_appearance()

    def _execute_drop(self, target_plot: 'SpectralPlotWindow') -> bool:
        """
        Execute the drop operation.

        Args:
            target_plot: The plot to drop the spectrum on

        Returns:
            True if drop was successful
        """
        if self._dragged_spectrum is None:
            return False

        return target_plot.accept_drop(self._dragged_spectrum)

    def cancel_drag(self) -> None:
        """Cancel the current drag operation and clean up."""
        self._is_dragging = False
        self._source_plot = None
        self._dragged_spectrum = None

        # Destroy the indicator
        if self._drag_indicator:
            try:
                self._drag_indicator.destroy()
            except tk.TclError:
                pass
            self._drag_indicator = None

        # Unbind global events from root
        try:
            self._root.unbind("<B1-Motion>")
            self._root.unbind("<ButtonRelease-1>")
        except tk.TclError:
            pass

        # Unbind from all plot windows and restore appearances
        plot_manager = self._plot_manager()
        for plot in plot_manager.get_all_plots().values():
            if plot.is_window_valid():
                try:
                    plot.window.unbind("<B1-Motion>")
                    plot.window.unbind("<ButtonRelease-1>")
                except tk.TclError:
                    pass
                plot._update_active_appearance()

    def is_dragging(self) -> bool:
        """Check if a drag operation is in progress."""
        return self._is_dragging

    def reset(self) -> None:
        """Reset the manager (primarily for testing)."""
        self.cancel_drag()


# Global singleton instance
_drag_drop_manager: Optional[DragDropManager] = None


def get_drag_drop_manager(root: tk.Tk = None) -> DragDropManager:
    """Get the global DragDropManager singleton."""
    global _drag_drop_manager
    if _drag_drop_manager is None:
        if root is None:
            raise ValueError("Root window required for first instantiation")
        _drag_drop_manager = DragDropManager(root)
    return _drag_drop_manager
