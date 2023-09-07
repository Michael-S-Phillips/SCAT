import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.widgets import Cursor, SpanSelector
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import spectral
import numpy as np

class HyperspectralAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Hyperspectral Image Analyzer")

        self.default_rgb_bands = [29, 19, 9]  # Default bands for RGB display
        self.default_stretch = 'linear'  # Default stretch for RGB display
        self.create_main_ui()
        self.create_menu()
        self.spectral_window = None 

    def create_main_ui(self):
        self.main_ui_frame = tk.Frame(self.root)
        self.main_ui_frame.grid(row=0, column=1, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
        self.root.grid_rowconfigure(0, weight=1)  # Allow row 0 to expand vertically with window resize
        self.root.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize

        # Create canvas and place it in the main_ui_frame
        self.create_canvas()

        # Create a sub-frame for the buttons
        self.button_frame = tk.Frame(self.main_ui_frame)
        self.button_frame.grid(row=0, column=1, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Allow row 0 to expand vertically with window resize
        self.main_ui_frame.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize

        # Create band selection options
        self.red_band_label = tk.Label(self.button_frame, text="Red Band:")
        self.red_band_label.grid(row=0, column=0, sticky=tk.W)

        self.red_band_var = tk.StringVar()
        self.red_band_menu = ttk.Combobox(self.button_frame, textvariable=self.red_band_var, state="readonly")
        self.red_band_menu.grid(row=0, column=1, sticky=tk.W)

        self.green_band_label = tk.Label(self.button_frame, text="Green Band:")
        self.green_band_label.grid(row=1, column=0, sticky=tk.W)

        self.green_band_var = tk.StringVar()
        self.green_band_menu = ttk.Combobox(self.button_frame, textvariable=self.green_band_var, state="readonly")
        self.green_band_menu.grid(row=1, column=1, sticky=tk.W)

        self.blue_band_label = tk.Label(self.button_frame, text="Blue Band:")
        self.blue_band_label.grid(row=2, column=0, sticky=tk.W)

        self.blue_band_var = tk.StringVar()
        self.blue_band_menu = ttk.Combobox(self.button_frame, textvariable=self.blue_band_var, state="readonly")
        self.blue_band_menu.grid(row=2, column=1, sticky=tk.W)

        self.apply_button = tk.Button(self.button_frame, text="Apply", command=self.update_rgb_display)
        self.apply_button.grid(row=3, column=0, columnspan=2, sticky=tk.W)

        # Create RGB stretch value options
        self.red_stretch_label = tk.Label(self.button_frame, text="Red Stretch:")
        self.red_stretch_label.grid(row=4, column=0, sticky=tk.W)

        self.red_min_stretch_var = tk.DoubleVar(value=0.0)
        self.red_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.red_min_stretch_var)
        self.red_min_stretch_entry.grid(row=5, column=0, sticky=tk.W)

        self.red_max_stretch_var = tk.DoubleVar(value=1.0)
        self.red_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.red_max_stretch_var)
        self.red_max_stretch_entry.grid(row=5, column=1, sticky=tk.W)

        self.green_stretch_label = tk.Label(self.button_frame, text="Green Stretch:")
        self.green_stretch_label.grid(row=6, column=0, sticky=tk.W)

        self.green_min_stretch_var = tk.DoubleVar(value=0.0)
        self.green_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.green_min_stretch_var)
        self.green_min_stretch_entry.grid(row=7, column=0, sticky=tk.W)

        self.green_max_stretch_var = tk.DoubleVar(value=1.0)
        self.green_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.green_max_stretch_var)
        self.green_max_stretch_entry.grid(row=7, column=1, sticky=tk.W)

        self.blue_stretch_label = tk.Label(self.button_frame, text="Blue Stretch:")
        self.blue_stretch_label.grid(row=8, column=0, sticky=tk.W)

        self.blue_min_stretch_var = tk.DoubleVar(value=0.0)
        self.blue_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.blue_min_stretch_var)
        self.blue_min_stretch_entry.grid(row=9, column=0, sticky=tk.W)

        self.blue_max_stretch_var = tk.DoubleVar(value=1.0)
        self.blue_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.blue_max_stretch_var)
        self.blue_max_stretch_entry.grid(row=9, column=1, sticky=tk.W)

        self.apply_stretch_button = tk.Button(self.button_frame, text="Apply Stretch", command=self.update_rgb_display)
        self.apply_stretch_button.grid(row=10, column=0, columnspan=2, sticky=tk.W)

        # Configure rows and columns to expand with resizing
        self.button_frame.grid_columnconfigure(0, weight=1)  # Allow column 0 to expand horizontally with window resize
        self.button_frame.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize

    def create_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load Data", command=self.load_data)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        hitogram_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Plot Hitograms", menu=hitogram_menu)
        hitogram_menu.add_command(label="Plot Histogram", command=self.plot_histograms)

    def create_canvas(self):
        self.figure, self.ax = plt.subplots(figsize=(4, 4))
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.main_ui_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0)  # Use grid instead of pack

        # # Configure columns to expand with resizing
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Adjust the column number if needed
        self.main_ui_frame.grid_columnconfigure(0, weight=1)

        # Bind the click event to the canvas
        self.canvas.mpl_connect('button_press_event', self.on_canvas_click)

    def load_data(self):
        file_path = filedialog.askopenfilename(filetypes=[("ENVI Files", "*.hdr")])
        if file_path:
            self.data = spectral.io.envi.open(file_path)
            self.populate_wavelength_menus()
            self.display_data(self.default_rgb_bands)

    def populate_wavelength_menus(self):
        if self.data is not None:
            wavelengths = self.data.bands.centers
            self.red_band_menu["values"] = wavelengths
            self.green_band_menu["values"] = wavelengths
            self.blue_band_menu["values"] = wavelengths
            self.red_band_menu.set(wavelengths[self.default_rgb_bands[0]])
            self.green_band_menu.set(wavelengths[self.default_rgb_bands[1]])
            self.blue_band_menu.set(wavelengths[self.default_rgb_bands[2]])

    def update_rgb_display(self):
        try:
            red_band = float(self.red_band_var.get())
            green_band = float(self.green_band_var.get())
            blue_band = float(self.blue_band_var.get())
            red_index = self.data.bands.centers.index(red_band)
            green_index = self.data.bands.centers.index(green_band)
            blue_index = self.data.bands.centers.index(blue_band)
            self.display_data([red_index, green_index, blue_index])
        except ValueError:
            messagebox.showerror("Error", "Invalid wavelength. Please enter valid wavelengths.")

    def display_data(self, band_indices):

        red_stretch = (self.red_min_stretch_var.get(), self.red_max_stretch_var.get()) if hasattr(self, "red_min_stretch_var") else None
        green_stretch = (self.green_min_stretch_var.get(), self.green_max_stretch_var.get()) if hasattr(self, "green_min_stretch_var") else None
        blue_stretch = (self.blue_min_stretch_var.get(), self.blue_max_stretch_var.get()) if hasattr(self, "blue_min_stretch_var") else None

        rgb_image = self.data[:, :, [band_indices[0], band_indices[1], band_indices[2]]]

        if red_stretch is not None:
            rgb_image[:, :, 0] = self.stretch_band(rgb_image[:, :, 0], red_stretch)
        if green_stretch is not None:
            rgb_image[:, :, 1] = self.stretch_band(rgb_image[:, :, 1], green_stretch)
        if blue_stretch is not None:
            rgb_image[:, :, 2] = self.stretch_band(rgb_image[:, :, 2], blue_stretch)

        self.ax.imshow(np.array(rgb_image))
        self.canvas.draw()

    def stretch_band(self, band, stretch_range):
        min_val, max_val = stretch_range
        stretched_band = (band - min_val) / (max_val - min_val)
        stretched_band = np.clip(stretched_band, 0, 1) 
        return stretched_band
    
    def plot_histograms(self):
        if hasattr(self, "data"):
            red_band = float(self.red_band_var.get())
            green_band = float(self.green_band_var.get())
            blue_band = float(self.blue_band_var.get())
            red_index = self.data.bands.centers.index(red_band)
            green_index = self.data.bands.centers.index(green_band)
            blue_index = self.data.bands.centers.index(blue_band)

            red_channel = np.where(self.data[:, :, red_index] > 1, np.nan, self.data[:, :, red_index])
            green_channel = np.where(self.data[:, :, green_index] > 1, np.nan, self.data[:, :, green_index])
            blue_channel = np.where(self.data[:, :, blue_index] > 1, np.nan, self.data[:, :, blue_index])

            # Calculate min and max values across all three channels
            min_value = min(np.nanmin(red_channel), np.nanmin(green_channel), np.nanmin(blue_channel))
            max_value = max(np.nanmax(red_channel), np.nanmax(green_channel), np.nanmax(blue_channel))

            if not hasattr(self, "hist_figure"):
                # Create a new histogram plot window if it doesn't exist
                self.hist_figure, self.hist_axes = plt.subplots(3, 1, figsize=(8, 6))
                self.hist_figure.show()

                # Create the subplots within the figure
                self.hist_axes[0].hist([], bins=256, color='red', alpha=0.7, range=(min_value, max_value))
                self.hist_axes[0].set_title(f'Red Channel Histogram: {red_band}')
                self.hist_axes[0].legend()

                self.hist_axes[1].hist([], bins=256, color='green', alpha=0.7, range=(min_value, max_value))
                self.hist_axes[1].set_title(f'Green Channel Histogram: {green_band}')
                self.hist_axes[1].legend()

                self.hist_axes[2].hist([], bins=256, color='blue', alpha=0.7, range=(min_value, max_value))
                self.hist_axes[2].set_title(f'Blue Channel Histogram: {blue_band}')
                self.hist_axes[2].legend()

            # Update the data in the existing subplots
            self.hist_axes[0].cla()
            self.hist_axes[0].hist(red_channel.ravel(), bins=256, color='red', alpha=0.7, range=(min_value, max_value))
            self.hist_axes[0].axvline(self.red_min_stretch_var.get(), color='red', linestyle='--', label='Min Stretch')
            self.hist_axes[0].axvline(self.red_max_stretch_var.get(), color='red', linestyle='--', label='Max Stretch')
            self.hist_axes[0].set_title(f'Red Channel Histogram: {red_band}')
            # self.hist_axes[0].legend()

            self.hist_axes[1].cla()
            self.hist_axes[1].hist(green_channel.ravel(), bins=256, color='green', alpha=0.7, range=(min_value, max_value))
            self.hist_axes[1].axvline(self.green_min_stretch_var.get(), color='green', linestyle='--', label='Min Stretch')
            self.hist_axes[1].axvline(self.green_max_stretch_var.get(), color='green', linestyle='--', label='Max Stretch')
            self.hist_axes[1].set_title(f'Green Channel Histogram: {green_band}')
            # self.hist_axes[1].legend()

            self.hist_axes[2].cla()
            self.hist_axes[2].hist(blue_channel.ravel(), bins=256, color='blue', alpha=0.7, range=(min_value, max_value))
            self.hist_axes[2].axvline(self.blue_min_stretch_var.get(), color='blue', linestyle='--', label='Min Stretch')
            self.hist_axes[2].axvline(self.blue_max_stretch_var.get(), color='blue', linestyle='--', label='Max Stretch')
            self.hist_axes[2].set_title(f'Blue Channel Histogram: {blue_band}')
            # self.hist_axes[2].legend()
            
            # Add spacing between subplots
            plt.subplots_adjust(hspace=0.4)

            # Update the existing figure
            self.hist_figure.canvas.draw()
            self.hist_figure.canvas.flush_events()

            # Create draggable span selectors
            self.create_histogram_span_selectors()

        else:
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first.")

    def update_red_stretch(self, val_min, val_max):
        self.red_min_stretch_var.set(val_min)
        self.red_max_stretch_var.set(val_max)
        self.plot_histograms()
        self.update_rgb_display()

    def update_green_stretch(self, val_min, val_max):
        self.green_min_stretch_var.set(val_min)
        self.green_max_stretch_var.set(val_max)
        self.plot_histograms()
        self.update_rgb_display()

    def update_blue_stretch(self, val_min, val_max):
        self.blue_min_stretch_var.set(val_min)
        self.blue_max_stretch_var.set(val_max)
        self.plot_histograms()
        self.update_rgb_display()

    def create_histogram_span_selectors(self):
        if not hasattr(self, "hist_span_selectors"):
            self.hist_span_selectors = []

        if not hasattr(self, "hist_figure"):
            return

        for selector in self.hist_span_selectors:
            selector.disconnect_events()

        self.hist_span_selectors.clear()

        red_selector = SpanSelector(self.hist_axes[0], self.update_red_stretch, 'horizontal', useblit=True)
        green_selector = SpanSelector(self.hist_axes[1], self.update_green_stretch, 'horizontal', useblit=True)
        blue_selector = SpanSelector(self.hist_axes[2], self.update_blue_stretch, 'horizontal', useblit=True)

        self.hist_span_selectors.extend([red_selector, green_selector, blue_selector])

    def on_canvas_click(self, event):
        if hasattr(self, "data"):
            if event.inaxes == self.ax:
                x, y = int(event.xdata), int(event.ydata)
                self.spectrum = self.data[y, x, :].flatten()
                self.spectrum = np.where(self.spectrum < 0, np.nan, self.spectrum)
                self.spectrum = np.where(self.spectrum > 1, np.nan, self.spectrum)
                self.update_spectral_plot()
        else:
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first.")


    def update_spectral_plot(self):
        if self.spectral_window is None or not self.spectral_window.winfo_exists():
            self.create_spectral_plot()
        else:
            self.spectral_line.set_ydata(self.spectrum)
            xmin, xmax = self.spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xmax))
            min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.spectral_canvas.draw()

    def create_spectral_plot(self):
        self.spectral_window = tk.Toplevel(self.root)
        self.spectral_window.title("Spectral Plot")
        
        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.spectral_window)
        ui_frame.pack(fill=tk.X)

        spectral_figure, self.spectral_ax = plt.subplots(figsize=(5,3))
        self.spectral_line, = self.spectral_ax.plot(self.data.bands.centers, self.spectrum)
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xmax))
        min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

        self.spectral_canvas = FigureCanvasTkAgg(spectral_figure, master=self.spectral_window)
        self.spectral_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset x-axis span
        self.reset_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_x_axis_span)
        self.reset_x_axis_button.pack(side=tk.RIGHT)

        # Add built-in span options to a dropdown menu
        span_options = ["Full Span", "410 - 1000 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
        self.span_var = tk.StringVar()
        self.span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(ui_frame, textvariable=self.span_var, values=span_options, state="readonly")
        span_menu.pack(side=tk.RIGHT)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_x_axis_span)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(self.spectral_canvas, self.spectral_window)
        toolbar.update()
        self.spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Add span selector for x-axis
        self.create_x_axis_span_selector(self.data.bands.centers)

    def reset_x_axis_span(self):
        # Reset x-axis span to the default range
        self.spectral_ax.set_xlim(self.data.bands.centers[0], self.data.bands.centers[-1])
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xmax))
        min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
        self.spectral_canvas.draw()

    def update_x_axis_span(self, event):
        selected_span = self.span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.data.bands.centers[0], self.data.bands.centers[-1]),
            "410 - 1000 nm": (410, 1000),
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900)
        }
        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            self.spectral_ax.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.data.bands.centers) - xlim[1]))

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
            self.spectral_ax, on_x_span_select, 'horizontal', useblit=True)


if __name__ == "__main__":
    root = tk.Tk()
    app = HyperspectralAnalyzer(root)
    root.mainloop()
