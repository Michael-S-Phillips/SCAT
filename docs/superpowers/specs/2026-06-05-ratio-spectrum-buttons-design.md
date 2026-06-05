# Dedicated Spectrum / Ratio Spectrum Plot Buttons

**Date:** 2026-06-05
**Status:** Approved

## Problem

The per-polygon "Plot Spectrum" and "Plot Ratio Spectrum" actions only exist in
the Treeview context menu, which is bound to `<Button-2>` (the X11 middle mouse
button). When SCAT runs on WSL and is accessed over a remote desktop from a
laptop, middle-click is awkward or impossible, so users cannot reach the
"Plot Ratio Spectrum" option.

There is an existing "Ratio Spectra" button, but it lives in the polygons
spectral plot window and plots *all* polygons' ratio spectra together — it is
not a substitute for plotting a single selected polygon's ratio spectrum.

## Goal

Provide an always-visible, mouse-friendly way to plot a single polygon's
spectrum and ratio spectrum, working identically to the existing context-menu
items.

## Design

### UI

Add two buttons to the Polygons Menu window's button bar (the same `ui_frame`
that already holds "Save ROIs", "Load ROIs", and the color dropdown, created in
`create_polygons_menu_window`):

- **`Plot Spectrum`** — calls `plot_selected_polygon(d=False)`
- **`Plot Ratio Spectrum`** — calls `plot_selected_polygon(d=True)`

The buttons act on the row currently selected in `self.polygon_table`. The user
already selects a row by clicking it (this highlights the polygon via
`highlight_polygon`), so no new selection mechanism is needed.

### Logic

Refactor `plot_selected_polygon(self, d=None)` to resolve the polygon index in
this order:

1. The current Treeview selection (`self.polygon_table.selection()`) — used by
   the new buttons.
2. Fallback to the existing `self.tmp_event` row identification — preserves the
   existing context-menu behavior.

If no row can be resolved (nothing selected and no usable event), show a
`messagebox.showwarning` telling the user to select a polygon row first, and
return without error.

Everything downstream is unchanged: the resolved `pc_index` is passed to
`spectral_window.create_spectral_plot(app, pc_index, d=d)`, so ratio/denominator
behavior is identical to today.

### Non-goals (explicitly out of scope)

- Keyboard shortcuts
- Adding/fixing the `<Button-3>` binding
- Modifier/double-click variants
- Removing or changing the existing context menu (it stays as-is)

## Testing

Manual verification (Tkinter GUI, no automated test harness in this project):

1. Load a cube, draw two polygons, mark one as denominator.
2. Select a non-denominator row, click **Plot Spectrum** → spectrum plot opens.
3. Select the same row, click **Plot Ratio Spectrum** → ratio spectrum plot
   opens (same result as the old context-menu item).
4. Click a button with no row selected → warning dialog, no crash.
5. Right-click (Button-2) context menu still works as before.
