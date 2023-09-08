import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.widgets import Cursor, SpanSelector
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import PickEvent
import spectral
import numpy as np

class HyperspectralAnalyzer:
    '''
    This class creates a GUI for analyzing hyperspectral data. It allows the user to load a hyperspectral and analyze an image 
    '''
    def __init__(self, root):
        '''
        Initializes the GUI and creates the main UI elements.
        '''
        self.root = root
        self.root.title("Hyperspectral Image Analyzer")

        self.default_rgb_bands = [29, 19, 9]  # Default bands for RGB display of spectral data
        self.default_parameter_bands = [0, 8, 9]  # Default bands for RGB display of parameter data
        self.default_stretch = 'linear'  # Default stretch for RGB display
        self.create_main_ui()
        self.create_menu()
        self.spectral_window = None 
        self.right_hist_window = None
        self.left_hist_window = None
        self.dragging = False
        self.drag_start_x = None
        self.drag_start_y = None

    def create_main_ui(self):
        '''
        Creates the main UI elements for the application.
        '''
        # Create a frame to hold UI elements with a fixed size
        self.main_ui_frame = tk.Frame(self.root)
        self.main_ui_frame.grid(row=0, column=1, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
        self.root.grid_rowconfigure(0, weight=1)  # Allow row 0 to expand vertically with window resize
        self.root.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize

        # Create canvases and place them in the main_ui_frame
        self.create_left_canvas()
        self.create_right_canvas()

        # Create a sub-frame for the left canvas buttons
        self.button_frame = tk.Frame(self.main_ui_frame)
        self.button_frame.grid(row=1, column=0, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Allow row 0 to expand vertically with window resize
        self.main_ui_frame.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize
        
        # Create band selection options for left canvas
        self.left_red_band_label = tk.Label(self.button_frame, text="Red Band:")
        self.left_red_band_label.grid(row=0, column=0, sticky=tk.W)

        self.left_red_band_var = tk.StringVar()
        self.left_red_band_menu = ttk.Combobox(self.button_frame, textvariable=self.left_red_band_var, state="readonly")
        self.left_red_band_menu.grid(row=0, column=1, sticky=tk.W)

        self.left_green_band_label = tk.Label(self.button_frame, text="Green Band:")
        self.left_green_band_label.grid(row=1, column=0, sticky=tk.W)

        self.left_green_band_var = tk.StringVar()
        self.left_green_band_menu = ttk.Combobox(self.button_frame, textvariable=self.left_green_band_var, state="readonly")
        self.left_green_band_menu.grid(row=1, column=1, sticky=tk.W)

        self.left_blue_band_label = tk.Label(self.button_frame, text="Blue Band:")
        self.left_blue_band_label.grid(row=2, column=0, sticky=tk.W)

        self.left_blue_band_var = tk.StringVar()
        self.left_blue_band_menu = ttk.Combobox(self.button_frame, textvariable=self.left_blue_band_var, state="readonly")
        self.left_blue_band_menu.grid(row=2, column=1, sticky=tk.W)

        self.apply_button = tk.Button(self.button_frame, text="Apply", command=self.update_left_display)
        self.apply_button.grid(row=3, column=0, columnspan=2, sticky=tk.W)

        # Create RGB stretch value options for left canvas
        self.left_red_stretch_label = tk.Label(self.button_frame, text="Red Stretch:")
        self.left_red_stretch_label.grid(row=4, column=0, sticky=tk.W)

        self.left_red_min_stretch_var = tk.DoubleVar(value=0.0)
        self.left_red_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.left_red_min_stretch_var)
        self.left_red_min_stretch_entry.grid(row=5, column=0, sticky=tk.W)

        self.left_red_max_stretch_var = tk.DoubleVar(value=1.0)
        self.left_red_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.left_red_max_stretch_var)
        self.left_red_max_stretch_entry.grid(row=5, column=1, sticky=tk.W)

        self.left_green_stretch_label = tk.Label(self.button_frame, text="Green Stretch:")
        self.left_green_stretch_label.grid(row=6, column=0, sticky=tk.W)

        self.left_green_min_stretch_var = tk.DoubleVar(value=0.0)
        self.left_green_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.left_green_min_stretch_var)
        self.left_green_min_stretch_entry.grid(row=7, column=0, sticky=tk.W)

        self.left_green_max_stretch_var = tk.DoubleVar(value=1.0)
        self.left_green_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.left_green_max_stretch_var)
        self.left_green_max_stretch_entry.grid(row=7, column=1, sticky=tk.W)

        self.left_blue_stretch_label = tk.Label(self.button_frame, text="Blue Stretch:")
        self.left_blue_stretch_label.grid(row=8, column=0, sticky=tk.W)

        self.left_blue_min_stretch_var = tk.DoubleVar(value=0.0)
        self.left_blue_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.left_blue_min_stretch_var)
        self.left_blue_min_stretch_entry.grid(row=9, column=0, sticky=tk.W)

        self.left_blue_max_stretch_var = tk.DoubleVar(value=1.0)
        self.left_blue_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.left_blue_max_stretch_var)
        self.left_blue_max_stretch_entry.grid(row=9, column=1, sticky=tk.W)

        self.apply_stretch_button = tk.Button(self.button_frame, text="Apply Stretch", command=self.update_left_display)
        self.apply_stretch_button.grid(row=10, column=0, columnspan=2, sticky=tk.W)

        # 
        # Create a sub-frame for the right canvas buttons
        self.button_frame = tk.Frame(self.main_ui_frame)
        self.button_frame.grid(row=1, column=1, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Allow row 0 to expand vertically with window resize
        self.main_ui_frame.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize

        # Create band selection options for the right canvas
        self.right_red_band_label = tk.Label(self.button_frame, text="Red Band:")
        self.right_red_band_label.grid(row=0, column=0, sticky=tk.W)

        self.right_red_band_var = tk.StringVar()
        self.right_red_band_menu = ttk.Combobox(self.button_frame, textvariable=self.right_red_band_var, state="readonly")
        self.right_red_band_menu.grid(row=0, column=1, sticky=tk.W)

        self.right_green_band_label = tk.Label(self.button_frame, text="Green Band:")
        self.right_green_band_label.grid(row=1, column=0, sticky=tk.W)

        self.right_green_band_var = tk.StringVar()
        self.right_green_band_menu = ttk.Combobox(self.button_frame, textvariable=self.right_green_band_var, state="readonly")
        self.right_green_band_menu.grid(row=1, column=1, sticky=tk.W)

        self.right_blue_band_label = tk.Label(self.button_frame, text="Blue Band:")
        self.right_blue_band_label.grid(row=2, column=0, sticky=tk.W)

        self.right_blue_band_var = tk.StringVar()
        self.right_blue_band_menu = ttk.Combobox(self.button_frame, textvariable=self.right_blue_band_var, state="readonly")
        self.right_blue_band_menu.grid(row=2, column=1, sticky=tk.W)

        self.apply_button = tk.Button(self.button_frame, text="Apply", command=self.update_right_display)
        self.apply_button.grid(row=3, column=0, columnspan=2, sticky=tk.W)

        # Create RGB stretch value options for the right canvas
        self.right_red_stretch_label = tk.Label(self.button_frame, text="Red Stretch:")
        self.right_red_stretch_label.grid(row=4, column=0, sticky=tk.W)

        self.right_red_min_stretch_var = tk.DoubleVar(value=0.0)
        self.right_red_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.right_red_min_stretch_var)
        self.right_red_min_stretch_entry.grid(row=5, column=0, sticky=tk.W)

        self.right_red_max_stretch_var = tk.DoubleVar(value=1.0)
        self.right_red_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.right_red_max_stretch_var)
        self.right_red_max_stretch_entry.grid(row=5, column=1, sticky=tk.W)

        self.right_green_stretch_label = tk.Label(self.button_frame, text="Green Stretch:")
        self.right_green_stretch_label.grid(row=6, column=0, sticky=tk.W)

        self.right_green_min_stretch_var = tk.DoubleVar(value=0.0)
        self.right_green_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.right_green_min_stretch_var)
        self.right_green_min_stretch_entry.grid(row=7, column=0, sticky=tk.W)

        self.right_green_max_stretch_var = tk.DoubleVar(value=1.0)
        self.right_green_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.right_green_max_stretch_var)
        self.right_green_max_stretch_entry.grid(row=7, column=1, sticky=tk.W)

        self.right_blue_stretch_label = tk.Label(self.button_frame, text="Blue Stretch:")
        self.right_blue_stretch_label.grid(row=8, column=0, sticky=tk.W)

        self.right_blue_min_stretch_var = tk.DoubleVar(value=0.0)
        self.right_blue_min_stretch_entry = tk.Entry(self.button_frame, textvariable=self.right_blue_min_stretch_var)
        self.right_blue_min_stretch_entry.grid(row=9, column=0, sticky=tk.W)

        self.right_blue_max_stretch_var = tk.DoubleVar(value=1.0)
        self.right_blue_max_stretch_entry = tk.Entry(self.button_frame, textvariable=self.right_blue_max_stretch_var)
        self.right_blue_max_stretch_entry.grid(row=9, column=1, sticky=tk.W)

        self.apply_stretch_button = tk.Button(self.button_frame, text="Apply Stretch", command=self.update_right_display)
        self.apply_stretch_button.grid(row=10, column=0, columnspan=2, sticky=tk.W)

        # Configure rows and columns to expand with resizing
        self.button_frame.grid_columnconfigure(0, weight=1)  # Allow column 0 to expand horizontally with window resize
        self.button_frame.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize

    def create_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load Hyperspectral Cube", command=self.load_left_data)
        file_menu.add_command(label="Load Band Parameter Image", command=self.load_right_data)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        hitogram_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Plot Hitograms", menu=hitogram_menu)
        hitogram_menu.add_command(label="Plot Left Frame Histogram", command=self.plot_left_histograms)
        hitogram_menu.add_command(label="Plot Right Frame Histogram", command=self.plot_right_histograms)

    def create_left_canvas(self):
        self.left_figure, self.left_ax = plt.subplots(figsize=(3, 3))
        self.left_canvas = FigureCanvasTkAgg(self.left_figure, master=self.main_ui_frame)
        self.left_canvas.get_tk_widget().grid(row=0, column=0)  # Use grid instead of pack

        # # Configure columns to expand with resizing
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Adjust the column number if needed
        self.main_ui_frame.grid_columnconfigure(0, weight=1)

        # Bind the click event to the canvas
        self.left_canvas.mpl_connect('button_press_event', self.on_left_canvas_click)
        self.left_canvas.mpl_connect('scroll_event', self.on_scroll)
        self.left_canvas.mpl_connect('button_release_event', self.on_canvas_release)
        self.left_canvas.mpl_connect('motion_notify_event', self.on_canvas_motion)

    def create_right_canvas(self):
        self.right_figure, self.right_ax = plt.subplots(figsize=(3, 3))
        self.right_canvas = FigureCanvasTkAgg(self.right_figure, master=self.main_ui_frame)
        self.right_canvas.get_tk_widget().grid(row=0, column=1)

        # Configure columns to expand with resizing
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Adjust the column number if needed
        self.main_ui_frame.grid_columnconfigure(1, weight=1)

        # Bind the click event to the canvas
        self.right_canvas.mpl_connect('button_press_event', self.on_left_canvas_click)
        self.right_canvas.mpl_connect('scroll_event', self.on_scroll)
        self.right_canvas.mpl_connect('button_release_event', self.on_canvas_release)
        self.right_canvas.mpl_connect('motion_notify_event', self.on_canvas_motion)

    def load_left_data(self):
        file_path = filedialog.askopenfilename(filetypes=[("ENVI Files", "*.hdr")])
        if file_path:
            self.left_data = spectral.io.envi.open(file_path)
            self.populate_left_wavelength_menus()
            self.display_left_data(self.default_rgb_bands)

    def load_right_data(self):
        file_path = filedialog.askopenfilename(filetypes=[("ENVI Files", "*.hdr")])
        if file_path:
            self.right_data = spectral.io.envi.open(file_path)
            self.populate_right_wavelength_menus()
            self.display_right_data(self.default_parameter_bands)

    def populate_left_wavelength_menus(self):
        if self.left_data is not None:
            wavelengths = self.left_data.bands.centers
            self.left_red_band_menu["values"] = wavelengths
            self.left_green_band_menu["values"] = wavelengths
            self.left_blue_band_menu["values"] = wavelengths
            self.left_red_band_menu.set(wavelengths[self.default_rgb_bands[0]])
            self.left_green_band_menu.set(wavelengths[self.default_rgb_bands[1]])
            self.left_blue_band_menu.set(wavelengths[self.default_rgb_bands[2]])

    def populate_right_wavelength_menus(self):
        if self.right_data is not None:
            wavelengths = self.right_data.metadata['band names']
            self.right_red_band_menu["values"] = wavelengths
            self.right_green_band_menu["values"] = wavelengths
            self.right_blue_band_menu["values"] = wavelengths
            self.right_red_band_menu.set(wavelengths[self.default_parameter_bands[0]])
            self.right_green_band_menu.set(wavelengths[self.default_parameter_bands[1]])
            self.right_blue_band_menu.set(wavelengths[self.default_parameter_bands[2]])

    def update_left_display(self):
        try:
            red_band = float(self.left_red_band_var.get())
            green_band = float(self.left_green_band_var.get())
            blue_band = float(self.left_blue_band_var.get())
            red_index = self.left_data.bands.centers.index(red_band)
            green_index = self.left_data.bands.centers.index(green_band)
            blue_index = self.left_data.bands.centers.index(blue_band)
            self.left_band_indices = [red_index, green_index, blue_index]
            self.display_left_data(self.left_band_indices)
        except ValueError:
            messagebox.showerror("Error", "Invalid wavelength. Please enter valid wavelengths.")

    def update_right_display(self):
        try:
            red_band = self.right_red_band_var.get()
            green_band = self.right_green_band_var.get()
            blue_band = self.right_blue_band_var.get()
            red_index = self.right_data.metadata['band names'].index(red_band)
            green_index = self.right_data.metadata['band names'].index(green_band)
            blue_index = self.right_data.metadata['band names'].index(blue_band)
            self.right_band_indices = [red_index, green_index, blue_index]
            self.display_right_data(self.right_band_indices)
        except ValueError:
            messagebox.showerror("Error", "Invalid wavelength. Please enter valid wavelengths.")

    def display_left_data(self, band_indices):

        left_red_stretch = (self.left_red_min_stretch_var.get(), self.left_red_max_stretch_var.get()) if hasattr(self, "left_red_min_stretch_var") else None
        left_green_stretch = (self.left_green_min_stretch_var.get(), self.left_green_max_stretch_var.get()) if hasattr(self, "left_green_min_stretch_var") else None
        left_blue_stretch = (self.left_blue_min_stretch_var.get(), self.left_blue_max_stretch_var.get()) if hasattr(self, "left_blue_min_stretch_var") else None

        left_rgb_image = self.left_data[:, :, [band_indices[0], band_indices[1], band_indices[2]]]

        if left_red_stretch is not None:
            left_rgb_image[:, :, 0] = self.stretch_band(left_rgb_image[:, :, 0], left_red_stretch)
        if left_green_stretch is not None:
            left_rgb_image[:, :, 1] = self.stretch_band(left_rgb_image[:, :, 1], left_green_stretch)
        if left_blue_stretch is not None:
            left_rgb_image[:, :, 2] = self.stretch_band(left_rgb_image[:, :, 2], left_blue_stretch)

        self.left_ax.imshow(np.array(left_rgb_image))
        self.left_canvas.draw()

    def display_right_data(self, band_indices):

        right_red_stretch = (self.right_red_min_stretch_var.get(), self.right_red_max_stretch_var.get()) if hasattr(self, "right_red_min_stretch_var") else None
        right_green_stretch = (self.right_green_min_stretch_var.get(), self.right_green_max_stretch_var.get()) if hasattr(self, "right_green_min_stretch_var") else None
        right_blue_stretch = (self.right_blue_min_stretch_var.get(), self.right_blue_max_stretch_var.get()) if hasattr(self, "right_blue_min_stretch_var") else None

        right_rgb_image = self.right_data[:, :, [band_indices[0], band_indices[1], band_indices[2]]]

        if right_red_stretch is not None:
            right_rgb_image[:, :, 0] = self.stretch_band(right_rgb_image[:, :, 0], right_red_stretch)
        if right_green_stretch is not None:
            right_rgb_image[:, :, 1] = self.stretch_band(right_rgb_image[:, :, 1], right_green_stretch)
        if right_blue_stretch is not None:
            right_rgb_image[:, :, 2] = self.stretch_band(right_rgb_image[:, :, 2], right_blue_stretch)

        self.right_ax.imshow(np.array(right_rgb_image))
        self.right_canvas.draw()

    def stretch_band(self, band, stretch_range):
        min_val, max_val = stretch_range
        stretched_band = (band - min_val) / (max_val - min_val)
        stretched_band = np.clip(stretched_band, 0, 1) 
        return stretched_band

    def plot_left_histograms(self):
        if hasattr(self, "left_data"):
            left_red_band = float(self.left_red_band_var.get())
            left_green_band = float(self.left_green_band_var.get())
            left_blue_band = float(self.left_blue_band_var.get())
            left_red_index = self.left_data.bands.centers.index(left_red_band)
            left_green_index = self.left_data.bands.centers.index(left_green_band)
            left_blue_index = self.left_data.bands.centers.index(left_blue_band)

            left_red_channel = np.where(self.left_data[:, :, left_red_index] > 1, np.nan, self.left_data[:, :, left_red_index])
            left_green_channel = np.where(self.left_data[:, :, left_green_index] > 1, np.nan, self.left_data[:, :, left_green_index])
            left_blue_channel = np.where(self.left_data[:, :, left_blue_index] > 1, np.nan, self.left_data[:, :, left_blue_index])

            # Calculate min and max values across all three channels
            left_min_value = min(np.nanmin(left_red_channel), np.nanmin(left_green_channel), np.nanmin(left_blue_channel))
            left_max_value = max(np.nanmax(left_red_channel), np.nanmax(left_green_channel), np.nanmax(left_blue_channel))

            if self.left_hist_window is None or not self.left_hist_window.winfo_exists():
                hist_args = left_min_value, left_max_value, left_red_band, left_green_band, left_blue_band
                self.create_left_histogram(hist_args)

            # Update the data in the existing subplots
            self.left_hist_axes[0].cla()
            self.left_hist_axes[0].hist(left_red_channel.ravel(), bins=256, color='red', alpha=0.7, range=(left_min_value, left_max_value))
            self.left_hist_axes[0].axvline(self.left_red_min_stretch_var.get(), color='red', linestyle='--', label='Min Stretch')
            self.left_hist_axes[0].axvline(self.left_red_max_stretch_var.get(), color='red', linestyle='--', label='Max Stretch')
            self.left_hist_axes[0].set_title(f'Left Red Channel Histogram: {left_red_band}')
            # self.left_hist_axes[0].legend()

            self.left_hist_axes[1].cla()
            self.left_hist_axes[1].hist(left_green_channel.ravel(), bins=256, color='green', alpha=0.7, range=(left_min_value, left_max_value))
            self.left_hist_axes[1].axvline(self.left_green_min_stretch_var.get(), color='green', linestyle='--', label='Min Stretch')
            self.left_hist_axes[1].axvline(self.left_green_max_stretch_var.get(), color='green', linestyle='--', label='Max Stretch')
            self.left_hist_axes[1].set_title(f'Left Green Channel Histogram: {left_green_band}')
            # self.left_hist_axes[1].legend()

            self.left_hist_axes[2].cla()
            self.left_hist_axes[2].hist(left_blue_channel.ravel(), bins=256, color='blue', alpha=0.7, range=(left_min_value, left_max_value))
            self.left_hist_axes[2].axvline(self.left_blue_min_stretch_var.get(), color='blue', linestyle='--', label='Min Stretch')
            self.left_hist_axes[2].axvline(self.left_blue_max_stretch_var.get(), color='blue', linestyle='--', label='Max Stretch')
            self.left_hist_axes[2].set_title(f'Left Blue Channel Histogram: {left_blue_band}')
            # self.left_hist_axes[2].legend()

            # Update the existing figure
            self.left_hist_figure.canvas.draw()
            self.left_hist_figure.canvas.flush_events()

            # Create draggable span selectors
            self.create_left_histogram_span_selectors()

        else:
            messagebox.showwarning("Warning", "No left frame data loaded. Load data into the left frame first.")

    def plot_right_histograms(self):
        if hasattr(self, "right_data"):
            right_red_band = self.right_red_band_var.get()
            right_green_band =self.right_green_band_var.get()
            right_blue_band = self.right_blue_band_var.get()
            right_red_index = self.right_data.metadata['band names'].index(right_red_band)
            right_green_index = self.right_data.metadata['band names'].index(right_green_band)
            right_blue_index = self.right_data.metadata['band names'].index(right_blue_band)

            right_red_channel = np.where(self.right_data[:, :, right_red_index] > 1, np.nan, self.right_data[:, :, right_red_index])
            right_green_channel = np.where(self.right_data[:, :, right_green_index] > 1, np.nan, self.right_data[:, :, right_green_index])
            right_blue_channel = np.where(self.right_data[:, :, right_blue_index] > 1, np.nan, self.right_data[:, :, right_blue_index])

            # Calculate min and max values across all three channels
            right_min_value = min(np.nanmin(right_red_channel), np.nanmin(right_green_channel), np.nanmin(right_blue_channel))
            right_max_value = max(np.nanmax(right_red_channel), np.nanmax(right_green_channel), np.nanmax(right_blue_channel))

            if self.right_hist_window is None or not self.right_hist_window.winfo_exists():
                hist_args = right_min_value, right_max_value, right_red_band, right_green_band, right_blue_band
                self.create_right_histogram(hist_args)

            # Update the data in the existing subplots
            self.right_hist_axes[0].cla()
            self.right_hist_axes[0].hist(right_red_channel.ravel(), bins=256, color='red', alpha=0.7, range=(right_min_value, right_max_value))
            self.right_hist_axes[0].axvline(self.right_red_min_stretch_var.get(), color='red', linestyle='--', label='Min Stretch')
            self.right_hist_axes[0].axvline(self.right_red_max_stretch_var.get(), color='red', linestyle='--', label='Max Stretch')
            self.right_hist_axes[0].set_title(f'Right Red Channel Histogram: {right_red_band}')
            # self.right_hist_axes[0].legend()

            self.right_hist_axes[1].cla()
            self.right_hist_axes[1].hist(right_green_channel.ravel(), bins=256, color='green', alpha=0.7, range=(right_min_value, right_max_value))
            self.right_hist_axes[1].axvline(self.right_green_min_stretch_var.get(), color='green', linestyle='--', label='Min Stretch')
            self.right_hist_axes[1].axvline(self.right_green_max_stretch_var.get(), color='green', linestyle='--', label='Max Stretch')
            self.right_hist_axes[1].set_title(f'Right Green Channel Histogram: {right_green_band}')
            # self.right_hist_axes[1].legend()

            self.right_hist_axes[2].cla()
            self.right_hist_axes[2].hist(right_blue_channel.ravel(), bins=256, color='blue', alpha=0.7, range=(right_min_value, right_max_value))
            self.right_hist_axes[2].axvline(self.right_blue_min_stretch_var.get(), color='blue', linestyle='--', label='Min Stretch')
            self.right_hist_axes[2].axvline(self.right_blue_max_stretch_var.get(), color='blue', linestyle='--', label='Max Stretch')
            self.right_hist_axes[2].set_title(f'Right Blue Channel Histogram: {right_blue_band}')
            # self.right_hist_axes[2].legend()

            # Update the existing figure
            self.right_hist_figure.canvas.draw()
            self.right_hist_figure.canvas.flush_events()

            # Create draggable span selectors
            self.create_right_histogram_span_selectors()

        else:
            messagebox.showwarning("Warning", "No right frame data loaded. Load data into the right frame first.")

    def create_left_histogram(self, hist_args):
        left_min_value, left_max_value, left_red_band, left_green_band, left_blue_band = hist_args
        self.left_hist_window = tk.Toplevel(self.root)
        self.left_hist_window.title("Hyperspectral Image Histogram")
        self.left_hist_window.geometry("600x800")
        
        # Create a frame to hold UI elements with a fixed size
        left_hist_ui_frame = tk.Frame(self.left_hist_window)
        left_hist_ui_frame.pack(fill=tk.X)

        # Create a new histogram plot window if it doesn't exist
        self.left_hist_figure, self.left_hist_axes = plt.subplots(3, 1, figsize=(4, 6))

        self.left_hist_canvas = FigureCanvasTkAgg(self.left_hist_figure, master=self.left_hist_window)
        self.left_hist_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset the histogram
        self.reset_x_axis_button = tk.Button(left_hist_ui_frame, text="Reset Histogram", command=self.plot_left_histograms)
        self.reset_x_axis_button.pack(side=tk.RIGHT)

        # Create the subplots within the figure
        self.left_hist_axes[0].hist([], bins=256, color='red', alpha=0.7, range=(left_min_value, left_max_value))
        self.left_hist_axes[0].set_title(f'Left Red Channel Histogram: {left_red_band}')
        self.left_hist_axes[0].legend()

        self.left_hist_axes[1].hist([], bins=256, color='green', alpha=0.7, range=(left_min_value, left_max_value))
        self.left_hist_axes[1].set_title(f'Left Green Channel Histogram: {left_green_band}')
        self.left_hist_axes[1].legend()

        self.left_hist_axes[2].hist([], bins=256, color='blue', alpha=0.7, range=(left_min_value, left_max_value))
        self.left_hist_axes[2].set_title(f'Left Blue Channel Histogram: {left_blue_band}')
        self.left_hist_axes[2].legend()
            
        # Add spacing between subplots
        self.left_hist_figure.subplots_adjust(hspace=0.8)

        # Create a toolbar for the spectral plot
        left_hist_toolbar = NavigationToolbar2Tk(self.left_hist_canvas, self.left_hist_window)
        left_hist_toolbar.update()
        self.left_hist_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def create_right_histogram(self, hist_args):
        right_min_value, right_max_value, right_red_band, right_green_band, right_blue_band = hist_args
        self.right_hist_window = tk.Toplevel(self.root)
        self.right_hist_window.title("Band Parameter Histogram")
        # self.right_hist_window.geometry("600x800")

        # Create a frame to hold UI elements with a fixed size
        right_hist_ui_frame = tk.Frame(self.right_hist_window)
        right_hist_ui_frame.pack(fill=tk.X)

        # Create a new histogram plot window if it doesn't exist
        self.right_hist_figure, self.right_hist_axes = plt.subplots(3, 1, figsize=(5, 6))

        self.right_hist_canvas = FigureCanvasTkAgg(self.right_hist_figure, master=self.right_hist_window)
        self.right_hist_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset the histogram
        self.reset_x_axis_button = tk.Button(right_hist_ui_frame, text="Reset Histogram", command=self.plot_right_histograms)
        self.reset_x_axis_button.pack(side=tk.RIGHT)

        # Create the subplots within the figure
        self.right_hist_axes[0].hist([], bins=256, color='red', alpha=0.7, range=(right_min_value, right_max_value))
        self.right_hist_axes[0].set_title(f'Right Red Channel Histogram: {right_red_band}', fontsize=4)
        self.right_hist_axes[0].legend()

        self.right_hist_axes[1].hist([], bins=256, color='green', alpha=0.7, range=(right_min_value, right_max_value))
        self.right_hist_axes[1].set_title(f'Right Green Channel Histogram: {right_green_band}', fontsize=4)
        self.right_hist_axes[1].legend()

        self.right_hist_axes[2].hist([], bins=256, color='blue', alpha=0.7, range=(right_min_value, right_max_value))
        self.right_hist_axes[2].set_title(f'Right Blue Channel Histogram: {right_blue_band}', fontsize=4)
        self.right_hist_axes[2].legend()
            
        # Add spacing between subplots
        self.right_hist_figure.subplots_adjust(hspace=0.8)

        # Create a toolbar for the spectral plot
        right_hist_toolbar = NavigationToolbar2Tk(self.right_hist_canvas, self.right_hist_window)
        right_hist_toolbar.update()
        self.right_hist_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

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

        left_red_selector = SpanSelector(self.left_hist_axes[0], self.update_left_red_stretch, 'horizontal', useblit=True)
        left_green_selector = SpanSelector(self.left_hist_axes[1], self.update_left_green_stretch, 'horizontal', useblit=True)
        left_blue_selector = SpanSelector(self.left_hist_axes[2], self.update_left_blue_stretch, 'horizontal', useblit=True)

        self.left_hist_span_selectors.extend([left_red_selector, left_green_selector, left_blue_selector])

    def create_right_histogram_span_selectors(self):
        if not hasattr(self, "right_hist_span_selectors"):
            self.right_hist_span_selectors = []

        if not hasattr(self, "right_hist_figure"):
            return

        for selector in self.right_hist_span_selectors:
            selector.disconnect_events()

        self.right_hist_span_selectors.clear()

        right_red_selector = SpanSelector(self.right_hist_axes[0], self.update_right_red_stretch, 'horizontal', useblit=True)
        right_green_selector = SpanSelector(self.right_hist_axes[1], self.update_right_green_stretch, 'horizontal', useblit=True)
        right_blue_selector = SpanSelector(self.right_hist_axes[2], self.update_right_blue_stretch, 'horizontal', useblit=True)

        self.right_hist_span_selectors.extend([right_red_selector, right_green_selector, right_blue_selector])

    def on_left_canvas_click(self, event):
        if hasattr(self, "left_data"):
            if event.inaxes == self.left_ax or event.inaxes == self.right_ax:
                if event.button == 1: #left mouse button press
                    self.drag_start_x = event.xdata
                    self.drag_start_y = event.ydata
                    if not self.dragging:
                        # single click for plotting spectra
                        x, y = int(event.xdata), int(event.ydata)
                        self.spectrum = self.left_data[y, x, :].flatten()
                        self.spectrum = np.where(self.spectrum < 0, np.nan, self.spectrum)
                        self.spectrum = np.where(self.spectrum > 1, np.nan, self.spectrum)
                        self.update_spectral_plot()
                    self.dragging = True
        elif event.button == 3:  # Right mouse button press
            # Implement additional functionality if needed
            pass
        else:
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first into left frame.")

    def on_canvas_release(self, event):
        if self.dragging:
            self.dragging = False
            self.left_canvas.draw()
            self.right_canvas.draw()

    def on_canvas_motion(self, event):
        if self.dragging:
            if self.drag_start_x is not None and self.drag_start_y is not None:
                dx = event.xdata - self.drag_start_x
                dy = event.ydata - self.drag_start_y
                self.drag_start_x = event.xdata
                self.drag_start_y = event.ydata

                # Adjust the axis limits to pan the image
                xlim = self.left_ax.get_xlim()
                ylim = self.left_ax.get_ylim()
                self.left_ax.set_xlim(xlim[0] - dx, xlim[1] - dx)
                self.left_ax.set_ylim(ylim[0] - dy, ylim[1] - dy)
                self.right_ax.set_xlim(xlim[0] - dx, xlim[1] - dx)
                self.right_ax.set_ylim(ylim[0] - dy, ylim[1] - dy)

                self.left_canvas.draw()
                self.right_canvas.draw()
                
    def on_scroll(self, event):
        if hasattr(self, "left_data") and not hasattr(self, "right_data"):
            if event.inaxes == self.left_ax:
                x_data, y_data = event.xdata, event.ydata
                x_lim, y_lim = self.left_ax.get_xlim(), self.left_ax.get_ylim()
                
                # Define the zoom factor (adjust as needed)
                zoom_factor = 1.1 if event.button == 'down' else 1 / 1.1  # Zoom in or out

                # Adjust the axis limits centered on the mouse cursor position
                new_x_lim = [x_data - (x_data - x_lim[0]) * zoom_factor, x_data + (x_lim[1] - x_data) * zoom_factor]
                new_y_lim = [y_data - (y_data - y_lim[0]) * zoom_factor, y_data + (y_lim[1] - y_data) * zoom_factor]
                
                self.left_ax.set_xlim(new_x_lim)
                self.left_ax.set_ylim(new_y_lim)
                
                event.canvas.draw()
        if hasattr(self, "right_data") and not hasattr(self, "left_data"):
            if event.inaxes == self.right_ax:
                x_data, y_data = event.xdata, event.ydata
                x_lim, y_lim = self.right_ax.get_xlim(), self.right_ax.get_ylim()
                
                # Define the zoom factor (adjust as needed)
                zoom_factor = 1.1 if event.button == 'down' else 1 / 1.1  # Zoom in or out

                # Adjust the axis limits centered on the mouse cursor position
                new_x_lim = [x_data - (x_data - x_lim[0]) * zoom_factor, x_data + (x_lim[1] - x_data) * zoom_factor]
                new_y_lim = [y_data - (y_data - y_lim[0]) * zoom_factor, y_data + (y_lim[1] - y_data) * zoom_factor]
                
                self.right_ax.set_xlim(new_x_lim)
                self.right_ax.set_ylim(new_y_lim)
                
                event.canvas.draw()
        if hasattr(self, "right_data") and hasattr(self, "left_data"):
            if event.inaxes == self.right_ax:
                # left axis adjust
                x_data, y_data = event.xdata, event.ydata
                x_lim, y_lim = self.left_ax.get_xlim(), self.left_ax.get_ylim()
                
                # Define the zoom factor (adjust as needed)
                zoom_factor = 1.1 if event.button == 'down' else 1 / 1.1  # Zoom in or out

                # Adjust the axis limits centered on the mouse cursor position
                new_x_lim = [x_data - (x_data - x_lim[0]) * zoom_factor, x_data + (x_lim[1] - x_data) * zoom_factor]
                new_y_lim = [y_data - (y_data - y_lim[0]) * zoom_factor, y_data + (y_lim[1] - y_data) * zoom_factor]
                
                self.left_ax.set_xlim(new_x_lim)
                self.left_ax.set_ylim(new_y_lim)
                self.right_ax.set_xlim(new_x_lim)
                self.right_ax.set_ylim(new_y_lim)
                
                event.canvas.draw()
                self.left_canvas.draw()
            elif event.inaxes == self.left_ax:
                # left axis adjust
                x_data, y_data = event.xdata, event.ydata
                x_lim, y_lim = self.left_ax.get_xlim(), self.left_ax.get_ylim()
                
                # Define the zoom factor (adjust as needed)
                zoom_factor = 1.1 if event.button == 'down' else 1 / 1.1  # Zoom in or out

                # Adjust the axis limits centered on the mouse cursor position
                new_x_lim = [x_data - (x_data - x_lim[0]) * zoom_factor, x_data + (x_lim[1] - x_data) * zoom_factor]
                new_y_lim = [y_data - (y_data - y_lim[0]) * zoom_factor, y_data + (y_lim[1] - y_data) * zoom_factor]
                
                self.left_ax.set_xlim(new_x_lim)
                self.left_ax.set_ylim(new_y_lim)

                # right axis adjust
                x_lim, y_lim = self.right_ax.get_xlim(), self.right_ax.get_ylim()
                
                # Define the zoom factor (adjust as needed)
                zoom_factor = 1.1 if event.button == 'down' else 1 / 1.1  # Zoom in or out

                # Adjust the axis limits centered on the mouse cursor position
                new_x_lim = [x_data - (x_data - x_lim[0]) * zoom_factor, x_data + (x_lim[1] - x_data) * zoom_factor]
                new_y_lim = [y_data - (y_data - y_lim[0]) * zoom_factor, y_data + (y_lim[1] - y_data) * zoom_factor]
                
                self.right_ax.set_xlim(new_x_lim)
                self.right_ax.set_ylim(new_y_lim)
                
                event.canvas.draw()
                self.right_canvas.draw()

    def update_spectral_plot(self):
        if self.spectral_window is None or not self.spectral_window.winfo_exists():
            self.create_spectral_plot()
        else:
            self.spectral_line.set_ydata(self.spectrum)
            xmin, xmax = self.spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xmax))
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
        self.spectral_line, = self.spectral_ax.plot(self.left_data.bands.centers, self.spectrum)
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xmax))
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
        self.create_x_axis_span_selector(self.left_data.bands.centers)

    def reset_x_axis_span(self):
        # Reset x-axis span to the default range
        self.spectral_ax.set_xlim(self.left_data.bands.centers[0], self.left_data.bands.centers[-1])
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xmax))
        min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
        self.spectral_canvas.draw()

    def update_x_axis_span(self, event):
        selected_span = self.span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_data.bands.centers[0], self.left_data.bands.centers[-1]),
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
            xmin_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_data.bands.centers) - xlim[1]))

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
