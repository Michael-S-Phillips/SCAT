# Dedicated Spectrum / Ratio Spectrum Plot Buttons Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add always-visible "Plot Spectrum" and "Plot Ratio Spectrum" buttons that act on the selected polygon table row, so ratio spectra can be plotted without the middle-click context menu (which is unusable over remote desktop).

**Architecture:** Refactor `plot_selected_polygon` to resolve the polygon index from the Treeview selection first (falling back to the existing mouse-event path), then add two buttons to the existing Polygons Menu button bar that call it with `d=False`/`d=True`.

**Tech Stack:** Python, Tkinter (`ttk.Treeview`, `tk.Button`), single-file app `scat.py`.

**Testing note:** This project has no automated test suite (monolithic Tkinter GUI). Verification is a Python import/syntax smoke check plus manual GUI steps from the spec.

---

### Task 1: Make `plot_selected_polygon` resolve index from table selection

**Files:**
- Modify: `scat.py` — `plot_selected_polygon` (currently around line 3391)

Current code:

```python
    def plot_selected_polygon(self, d=None):
        event = self.tmp_event
        item = self.polygon_table.identify_row(event.y)
        pc_index = int(self.polygon_table.item(item, "values")[0])
        ssw = spectral_window  # single_spectrum_window
        # ssw.create_window(app, pc_index) #old
        ssw.create_spectral_plot(app, pc_index, d=d)
        # spec = self.polygon_spectra[pc_index]
```

- [ ] **Step 1: Replace the method body to prefer the Treeview selection**

Replace the method above with:

```python
    def plot_selected_polygon(self, d=None):
        # Resolve the target row: prefer the current table selection (used by
        # the Plot Spectrum / Plot Ratio Spectrum buttons), then fall back to
        # the row under the last context-menu click (self.tmp_event).
        item = None
        selection = self.polygon_table.selection()
        if selection:
            item = selection[0]
        elif getattr(self, "tmp_event", None) is not None:
            item = self.polygon_table.identify_row(self.tmp_event.y)

        if not item:
            messagebox.showwarning(
                "No polygon selected",
                "Select a polygon row in the table first, then click the plot button.",
            )
            return

        pc_index = int(self.polygon_table.item(item, "values")[0])
        ssw = spectral_window  # single_spectrum_window
        ssw.create_spectral_plot(app, pc_index, d=d)
```

- [ ] **Step 2: Smoke-check that the file still imports / parses**

Run: `python -c "import ast; ast.parse(open('scat.py').read()); print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add scat.py
git commit -m "Resolve plot_selected_polygon target from table selection"
```

---

### Task 2: Add the two plot buttons to the Polygons Menu button bar

**Files:**
- Modify: `scat.py` — `create_polygons_menu_window`, just before `self.create_polygons_table()` (currently around line 2597-2601)

Current code (end of the button bar, in `ui_frame`):

```python
        self.polygon_color_dropdown.grid(
            row=self.polygons_menu_row, column=self.polygons_menu_col
        )
        self.polygons_menu_row += 1

        # ----------------------------------------------------------------
        # create a table to display polygon information
        self.create_polygons_table()
```

- [ ] **Step 1: Insert the two buttons between the dropdown row and the table**

Replace the block above with:

```python
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
```

- [ ] **Step 2: Smoke-check that the file still imports / parses**

Run: `python -c "import ast; ast.parse(open('scat.py').read()); print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add scat.py
git commit -m "Add Plot Spectrum / Plot Ratio Spectrum buttons to polygons menu"
```

---

### Task 3: Manual GUI verification

**Files:** none (manual)

- [ ] **Step 1: Launch and exercise the buttons**

Run: `python scat.py`

Verify (from the spec):
1. Load a cube, draw two polygons, mark one as denominator.
2. Open the Polygons Menu; the **Plot Spectrum** and **Plot Ratio Spectrum** buttons appear under the existing button bar.
3. Select a non-denominator row, click **Plot Spectrum** → spectrum plot opens.
4. Select the same row, click **Plot Ratio Spectrum** → ratio spectrum plot opens (same result as the old context-menu item).
5. Click a button with no row selected → warning dialog, no crash.
6. Right-click (Button-2) context menu still works as before.

(The author cannot run the GUI in this environment; the user performs this step.)
