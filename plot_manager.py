"""
PlotManager: Central registry for all spectral plot windows.

This module provides a singleton manager that tracks all plot windows,
handles active plot designation, and coordinates spectrum routing.
"""

import uuid
from typing import Optional, Dict, TYPE_CHECKING

if TYPE_CHECKING:
    from spectral_plot_window import SpectralPlotWindow


class PlotManager:
    """
    Central registry for managing spectral plot windows.

    Responsibilities:
    - Track all registered plot windows
    - Manage active plot designation
    - Route new spectra to the active plot
    - Clean up closed windows
    """

    _instance = None

    def __new__(cls):
        """Singleton pattern to ensure only one PlotManager exists."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if self._initialized:
            return
        self._initialized = True
        self._plots: Dict[str, 'SpectralPlotWindow'] = {}
        self._active_plot_id: Optional[str] = None
        self._on_active_changed_callbacks = []

    def register_plot(self, plot: 'SpectralPlotWindow') -> str:
        """
        Register a plot window and return its unique ID.

        Args:
            plot: The SpectralPlotWindow to register

        Returns:
            Unique string ID for the plot
        """
        plot_id = str(uuid.uuid4())[:8]
        self._plots[plot_id] = plot

        # If this is the first/only plot, make it active
        if len(self._plots) == 1 or self._active_plot_id is None:
            self.set_active_plot(plot_id)

        return plot_id

    def unregister_plot(self, plot_id: str) -> None:
        """
        Remove a plot from the registry (called on window close).

        Args:
            plot_id: The ID of the plot to unregister
        """
        if plot_id not in self._plots:
            return

        del self._plots[plot_id]

        # If the closed plot was active, activate another one
        if self._active_plot_id == plot_id:
            self._active_plot_id = None
            if self._plots:
                # Activate the most recently added plot
                next_plot_id = list(self._plots.keys())[-1]
                self.set_active_plot(next_plot_id)
            else:
                self._notify_active_changed(None)

    def set_active_plot(self, plot_id: str) -> bool:
        """
        Designate a plot as the active plot (receives new spectra).

        Args:
            plot_id: The ID of the plot to activate

        Returns:
            True if successful, False if plot_id not found
        """
        if plot_id not in self._plots:
            return False

        # Deactivate previous active plot
        if self._active_plot_id is not None and self._active_plot_id in self._plots:
            self._plots[self._active_plot_id].set_active(False)

        # Activate new plot
        self._active_plot_id = plot_id
        self._plots[plot_id].set_active(True)
        self._notify_active_changed(plot_id)

        return True

    def get_active_plot(self) -> Optional['SpectralPlotWindow']:
        """
        Get the currently active plot window.

        Returns:
            The active SpectralPlotWindow, or None if no plots exist
        """
        if self._active_plot_id is None:
            return None
        return self._plots.get(self._active_plot_id)

    def get_active_plot_id(self) -> Optional[str]:
        """Get the ID of the currently active plot."""
        return self._active_plot_id

    def get_plot(self, plot_id: str) -> Optional['SpectralPlotWindow']:
        """Get a plot by its ID."""
        return self._plots.get(plot_id)

    def get_all_plots(self) -> Dict[str, 'SpectralPlotWindow']:
        """Get all registered plots."""
        return self._plots.copy()

    def get_plot_count(self) -> int:
        """Get the number of registered plots."""
        return len(self._plots)

    def send_spectrum_to_active(self, spectrum_data: 'SpectrumData') -> bool:
        """
        Send a spectrum to the active plot.

        Args:
            spectrum_data: The SpectrumData to add to the active plot

        Returns:
            True if spectrum was sent, False if no active plot
        """
        active_plot = self.get_active_plot()
        if active_plot is None:
            return False

        active_plot.add_spectrum(spectrum_data)
        return True

    def on_active_changed(self, callback) -> None:
        """
        Register a callback to be called when the active plot changes.

        Args:
            callback: Function taking plot_id (or None) as argument
        """
        self._on_active_changed_callbacks.append(callback)

    def _notify_active_changed(self, plot_id: Optional[str]) -> None:
        """Notify all callbacks of active plot change."""
        for callback in self._on_active_changed_callbacks:
            try:
                callback(plot_id)
            except Exception:
                pass  # Don't let callback errors break the manager

    def cleanup_closed_windows(self) -> None:
        """Remove any plots whose windows have been closed."""
        closed_ids = []
        for plot_id, plot in self._plots.items():
            if not plot.is_window_valid():
                closed_ids.append(plot_id)

        for plot_id in closed_ids:
            self.unregister_plot(plot_id)

    def reset(self) -> None:
        """Reset the manager (primarily for testing)."""
        self._plots.clear()
        self._active_plot_id = None
        self._on_active_changed_callbacks.clear()


# Global singleton instance
_plot_manager: Optional[PlotManager] = None


def get_plot_manager() -> PlotManager:
    """Get the global PlotManager singleton."""
    global _plot_manager
    if _plot_manager is None:
        _plot_manager = PlotManager()
    return _plot_manager
