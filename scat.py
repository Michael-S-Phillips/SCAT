import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pickle
from matplotlib.widgets import SpanSelector
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import spectral
import numpy as np
import rasterio as rio
import cv2

import geopandas as gpd
from shapely.geometry import Polygon as sgp

class SpectralCubeAnalysisTool:
    '''
    TODO:
        - Adding more functionality to the polygon shape layer drawing:
            > naming polygons
            > better ability to edit individual polygons
            > ability to delete individual polygons
            > ability to cancel a polygon in the middle of drawing
            > a way to interact with polygons in a table format (edit names and colors)
            > speed up mean spectra plotting
        - Spectral plotting
            > add ability to save plots at figure quality
            > spectral smoothing routines (Rogers for CRISM, generic smoothing)
            > add ability to plot multiple spectra on the same plot for point and click spectra
            > add ability to plot library spectra
            > add ability to offset for clarity
        - Adding basic spectral analysis workflows, such as MNF and HyPyRameter
        - Add default Browse product combos for parameter display
        - Add ability to save statistics from parameter file over ROIs

    This class creates a GUI for analyzing hyperspectral data. It allows the user to load a hyperspectral and analyze an image 
    '''

    def __init__(self, root):
        '''
        Initializes the GUI and creates the main UI elements.
        '''
        self.root = root
        self.root.title("Hyperspectral Image Analyzer")

        self.default_rgb_bands = [29, 19, 9]  # Default bands for RGB display of spectral data
        self.default_parameter_bands = [15, 18, 19]  # Default bands for RGB display of parameter data, MAF
        self.default_stretch = 'linear'  # Default stretch for RGB display
        self.create_main_ui()
        self.create_menu()
        self.spectral_window = None 
        self.ratio_spectral_window = None 
        self.polygons_spectral_window = None
        self.right_hist_window = None
        self.left_hist_window = None
        self.right_is_parameter = False
        self.left_is_paramenter = False
        self.reset_right_hist = False
        self.reset_left_hist = False
        self.draw_polygons = False
        self.points = []
        self.current_polygon = []
        self.all_polygons = []
        self.polygon_colors = []
        self.polygon_spectra = []

    def create_main_ui(self):
        '''
        Creates the main UI elements for the application.
        '''
        # Create a frame to hold UI elements with a fixed size
        self.main_ui_frame = tk.Frame(self.root)
        self.main_ui_frame.grid(row=0, column=1, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
        self.root.grid_rowconfigure(0, weight=1)  # Allow row 0 to expand vertically with window resize
        self.root.grid_columnconfigure(1, weight=1)  # Allow column 1 to expand horizontally with window resize

        # Create a label to display coordinates
        self.coordinates_label = tk.Label(self.main_ui_frame, text="Lat: 0.0000, Lon: 0.0000")
        self.coordinates_label.grid(row=0, column=0, columnspan=2, sticky="nsew")
        # self.coordinates_label.grid_rowconfigure(row=0, weight=0)

        # frame for the buttons on the righthand side
        self.right_buttons_frame = tk.Frame(self.main_ui_frame)
        self.right_buttons_frame.grid(row=0, column=2, rowspan=2, sticky='new')
        self.main_ui_frame.grid_columnconfigure(2, weight=1)
        self.right_buttons_frame.grid_rowconfigure(0, weight=1)
        self.right_buttons_frame.grid_columnconfigure(0, weight=1)

        # Create a button to reset displays
        self.reset_display_button = tk.Button(self.right_buttons_frame, text="Reset Display Extent", command=self.reset_display_extent)
        self.reset_display_button.grid(row=0, column=0, sticky="new")
        # self.reset_display_button.grid_rowconfigure(row=0, weight=0)

        # Add a button to draw polygons
        self.draw_polygons_button = tk.Button(self.right_buttons_frame, text="Draw Polygons", command=self.create_polygons_menu_window)
        self.draw_polygons_button.grid(row=1, column=0, sticky = 'new')

        # Create a sub-frame for the left canvas buttons
        self.button_frame = tk.Frame(self.main_ui_frame)
        self.button_frame.grid(row=2, column=0, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
        self.main_ui_frame.grid_rowconfigure(0, weight=0)  # Allow row 0 to expand vertically with window resize
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

        # ----------------------------------------------------------------
        # Create a sub-frame for the right canvas buttons
        # ----------------------------------------------------------------
        self.button_frame = tk.Frame(self.main_ui_frame)
        self.button_frame.grid(row=2, column=1, sticky=tk.N + tk.S + tk.W + tk.E)  # Use grid and place it in column 1
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
        file_menu.add_command(label = "Load Session", command=self.load_state)
        file_menu.add_command(label = "Save", command=self.save_state)
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        hitogram_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Plot Hitograms", menu=hitogram_menu)
        hitogram_menu.add_command(label="Plot Left Frame Histogram", command=self.plot_left_histograms)
        hitogram_menu.add_command(label="Plot Right Frame Histogram", command=self.plot_right_histograms)

    def create_left_canvas(self):
        self.left_figure, self.left_ax = plt.subplots(figsize=(6, 6))
        self.left_ax.set_title("Hyperspectral Cube", fontsize= 8)
        self.left_ax.tick_params(axis='both', which='major', labelsize=8)
        # Set the margins and spacing to zero
        self.left_figure.tight_layout()
        self.left_figure.subplots_adjust(wspace=0, hspace=0)
        
        self.left_frame = tk.Frame(self.main_ui_frame)
        self.left_frame.grid(row=1, column=0, sticky="nsew")

        self.left_canvas = FigureCanvasTkAgg(self.left_figure, master=self.left_frame)
        self.left_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")  # Use grid instead of pack

        toolbar_frame = tk.Frame(self.left_frame)
        toolbar = NavigationToolbar2Tk(self.left_canvas, toolbar_frame)
        toolbar.update()
        toolbar_frame.grid(row=1, column=0, sticky="nsew")

        # Configure columns to expand with resizing
        self.left_frame.grid_rowconfigure(0, weight=1)  # Adjust the column number if needed
        self.left_frame.grid_columnconfigure(0, weight=1)
        self.left_canvas.get_tk_widget().grid_columnconfigure(0, weight=1)
        self.left_canvas.get_tk_widget().grid_rowconfigure(0, weight=1)

        if self.left_data.nrows/self.left_data.ncols > 2:
            # Create variables to store the current top-left pixel coordinates
            self.top_left_row = 0
            self.top_left_col = 0

            # Define the fixed size of the displayed portion
            self.left_display_rows = int(3*self.left_data.ncols)
            self.left_ax.set_ylim(int(3*self.left_data.ncols), 0)
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
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Adjust the column number if needed
        self.main_ui_frame.grid_rowconfigure(1, weight=1)  # Adjust the column number if needed
        self.main_ui_frame.grid_columnconfigure(0, weight=1)

        # Bind the click event to the canvas
        self.left_canvas.mpl_connect('button_press_event', self.on_left_canvas_click)
        self.left_canvas.mpl_connect('scroll_event', self.on_scroll)
        # self.left_canvas.mpl_connect('button_release_event', self.on_canvas_release)
        self.left_canvas.mpl_connect('motion_notify_event', self.on_canvas_motion)

    def create_right_canvas(self):
        self.right_figure, self.right_ax = plt.subplots(figsize=(6, 6))
        self.right_ax.set_title("Band Parameter Image", fontsize=8)
        self.right_ax.tick_params(axis='both', which='major', labelsize=8)
        # Set the margins and spacing to zero
        self.right_figure.tight_layout()
        self.right_figure.subplots_adjust(wspace=0, hspace=0)

        self.right_frame = tk.Frame(self.main_ui_frame)
        self.right_frame.grid(row=1, column=1, sticky="nsew")

        self.right_canvas = FigureCanvasTkAgg(self.right_figure, master=self.right_frame)
        self.right_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")  # Use grid instead of pack
        # self.right_canvas.get_tk_widget().grid(row=1, column=0)

        toolbar_frame = tk.Frame(self.right_frame)
        toolbar = NavigationToolbar2Tk(self.right_canvas, toolbar_frame)
        toolbar.update()
        toolbar_frame.grid(row=1, column=0, sticky="nsew")


        # Configure columns to expand with resizing
        self.right_frame.grid_rowconfigure(0, weight=1)  # Adjust the column number if needed
        self.right_frame.grid_columnconfigure(0, weight=1)
        self.right_canvas.get_tk_widget().grid_columnconfigure(0, weight=1)
        self.right_canvas.get_tk_widget().grid_rowconfigure(0, weight=1)

        if self.right_data.nrows/self.right_data.ncols > 2:
            # Create variables to store the current top-left pixel coordinates
            self.top_right_row = 0
            self.top_right_col = 0

            # Define the fixed size of the displayed portion
            self.right_display_rows = int(3*self.right_data.ncols)
            self.right_ax.set_ylim(int(3*self.right_data.ncols), 0)
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
        self.main_ui_frame.grid_rowconfigure(0, weight=1)  # Adjust the column number if needed
        self.main_ui_frame.grid_rowconfigure(1, weight=1)  # Adjust the column number if needed
        self.main_ui_frame.grid_columnconfigure(1, weight=1)

        # Bind the click event to the canvas
        self.right_canvas.mpl_connect('button_press_event', self.on_left_canvas_click)
        self.right_canvas.mpl_connect('scroll_event', self.on_scroll)
        # self.right_canvas.mpl_connect('button_release_event', self.on_canvas_release)
        self.right_canvas.mpl_connect('motion_notify_event', self.on_canvas_motion)

# ----------------------------------------------------------------
# loading the data
# ----------------------------------------------------------------
    def load_left_data(self, file_path = None):
        if not file_path:
            file_path = filedialog.askopenfilename(filetypes=[("ENVI Files", "*.hdr")])
        if file_path:
            self.left_data = spectral.io.envi.open(file_path)
            rio_path = file_path.replace("hdr", "img")
            self.left_rio = rio.open(rio_path)
            self.map_info = self.left_data.metadata['map info']
            self.crs = self.left_data.metadata['coordinate system string']
            self.x_dpp = 5/self.left_data.ncols
            self.y_dpp = 5/self.left_data.nrows
            if hasattr(self.left_data.metadata, 'projection info'):
                self.proj_info = self.left_data.metadata['projection info']
            else:
                self.proj_info = self.map_info[1:]
            self.left_transform = self.left_rio.transform
            self.transformer = rio.transform.AffineTransformer(self.left_transform)
            self.populate_left_wavelength_menus()
            self.create_left_canvas()
            self.create_left_scroll_bars()
            self.display_left_data(self.default_rgb_bands)

    def load_right_data(self, file_path = None):
        if not file_path:
            file_path = filedialog.askopenfilename(filetypes=[("ENVI Files", "*.hdr")])
        if file_path:
            self.right_data = spectral.io.envi.open(file_path)
            self.populate_right_wavelength_menus()
            self.create_right_canvas()
            self.create_right_scroll_bars()
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
                        if self.left_data.bands.centers[1] - self.left_data.bands.centers[-1] > 0:
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.left_data.bands.centers
                            self.left_wvl = sorted([float(i) for i in self.left_data.bands.centers])
                        else:
                            wavelengths = self.left_data.bands.centers
                            self.left_wvl = [float(i) for i in self.left_data.bands.centers]

                elif self.left_data.bands.centers is None and self.left_data.metadata['band names'] is not None:
                    if str(self.left_data.metadata['band names'][1]).isalpha():
                        # This is a band parameter image
                        self.left_is_parameter = True
                        wavelengths = self.left_data.metadata['band names']
                        self.left_wvl = self.left_data.metadata['band names']
                    else:
                        # This is a spectral cube
                        if self.left_data.metadata['band names'][1] - self.left_data.metadata['band names'][-1] > 0:
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.left_data.metadata['band names']
                            self.left_wvl = sorted([float(i) for i in self.left_data.metadata['band names']])
                        else:
                            # This is a band parameter image
                            wavelengths = self.left_data.metadata['band names']
                            self.left_wvl = self.left_data.metadata['band names']
            except:
                messagebox.showerror("Error", "Unable to load wavelength information.")

            self.left_red_band_menu["values"] = wavelengths
            self.left_green_band_menu["values"] = wavelengths
            self.left_blue_band_menu["values"] = wavelengths
            self.left_red_band_menu.set(wavelengths[self.default_rgb_bands[0]])
            self.left_green_band_menu.set(wavelengths[self.default_rgb_bands[1]])
            self.left_blue_band_menu.set(wavelengths[self.default_rgb_bands[2]])

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
                        if self.right_data.bands.centers[1] - self.right_data.bands.centers[-1] > 0:
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.right_data.bands.centers
                            self.right_wvl = sorted([float(i) for i in self.right_data.bands.centers])
                        else:
                            wavelengths = self.right_data.bands.centers
                            self.right_wvl = [float(i) for i in self.right_data.bands.centers]
                
                elif self.right_data.bands.centers is None and self.right_data.metadata['band names'] is not None:
                    if any(char.isalpha() for char in str(self.right_data.metadata['band names'][1])):
                        # this is a band parameter image
                        self.right_is_parameter = True
                        wavelengths = self.right_data.metadata['band names']
                        self.right_wvl = self.right_data.metadata['band names']
                    else:
                        # this is a spectral cube
                        if self.right_data.metadata['band names'][1] - self.right_data.metadata['band names'][-1] > 0:
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.right_data.metadata['band names']
                            self.right_wvl = sorted([float(i) for i in self.right_data.metadata['band names']])
                        else:
                            # this is a band parameter image
                            wavelengths = self.right_data.metadata['band names']
                            self.right_wvl = self.right_data.metadata['band names']
            except:
                messagebox.showerror("Error", "Unable to load wavelength information.")
            
            self.right_red_band_menu["values"] = wavelengths
            self.right_green_band_menu["values"] = wavelengths
            self.right_blue_band_menu["values"] = wavelengths
            self.right_red_band_menu.set(wavelengths[self.default_parameter_bands[0]])
            self.right_green_band_menu.set(wavelengths[self.default_parameter_bands[1]])
            self.right_blue_band_menu.set(wavelengths[self.default_parameter_bands[2]])

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
        except ValueError:
            messagebox.showerror("Error", "Invalid wavelength. Please enter valid wavelengths.")

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
        except ValueError:
            messagebox.showerror("Error", "Invalid wavelength. Please enter valid wavelengths.")

    def display_left_data(self, band_indices):

        left_red_stretch = (self.left_red_min_stretch_var.get(), self.left_red_max_stretch_var.get()) if hasattr(self, "left_red_min_stretch_var") else None
        left_green_stretch = (self.left_green_min_stretch_var.get(), self.left_green_max_stretch_var.get()) if hasattr(self, "left_green_min_stretch_var") else None
        left_blue_stretch = (self.left_blue_min_stretch_var.get(), self.left_blue_max_stretch_var.get()) if hasattr(self, "left_blue_min_stretch_var") else None

        # print(self.top_left_row, self.left_display_rows, self.top_left_col, self.left_display_cols)
        # left_rgb_image = self.left_data[self.top_left_row:self.top_left_row + self.left_display_rows,
        #                                 self.top_left_col:self.top_left_col + self.left_display_cols, 
        #                                 [band_indices[0], band_indices[1], band_indices[2]]]
        
        left_rgb_image = self.left_data[:,:,[band_indices[0], band_indices[1], band_indices[2]]]

        if left_red_stretch is not None:
            left_rgb_image[:, :, 0] = self.stretch_band(left_rgb_image[:, :, 0], left_red_stretch)
        if left_green_stretch is not None:
            left_rgb_image[:, :, 1] = self.stretch_band(left_rgb_image[:, :, 1], left_green_stretch)
        if left_blue_stretch is not None:
            left_rgb_image[:, :, 2] = self.stretch_band(left_rgb_image[:, :, 2], left_blue_stretch)

        self.left_ax.imshow(np.array(left_rgb_image))
        # show(np.transpose(left_rgb_image, (2,0,1)), ax=self.left_ax, transform=self.left_transform)
        self.left_canvas.draw()

    def display_right_data(self, band_indices):

        right_red_stretch = (self.right_red_min_stretch_var.get(), self.right_red_max_stretch_var.get()) if hasattr(self, "right_red_min_stretch_var") else None
        right_green_stretch = (self.right_green_min_stretch_var.get(), self.right_green_max_stretch_var.get()) if hasattr(self, "right_green_min_stretch_var") else None
        right_blue_stretch = (self.right_blue_min_stretch_var.get(), self.right_blue_max_stretch_var.get()) if hasattr(self, "right_blue_min_stretch_var") else None

        # right_rgb_image = self.right_data[self.top_right_row:self.top_right_row + self.right_display_rows,
        #                                 self.top_right_col:self.top_right_col + self.right_display_cols, 
        #                                 [band_indices[0], band_indices[1], band_indices[2]]]

        right_rgb_image = self.right_data[:,:,[band_indices[0], band_indices[1], band_indices[2]]]

        if right_red_stretch is not None:
            right_rgb_image[:, :, 0] = self.stretch_band(right_rgb_image[:, :, 0], right_red_stretch)
        if right_green_stretch is not None:
            right_rgb_image[:, :, 1] = self.stretch_band(right_rgb_image[:, :, 1], right_green_stretch)
        if right_blue_stretch is not None:
            right_rgb_image[:, :, 2] = self.stretch_band(right_rgb_image[:, :, 2], right_blue_stretch)

        self.right_ax.imshow(np.array(right_rgb_image))
        self.right_canvas.draw()
    
    def reset_display_extent(self):
        self.left_ax.set_ylim(self.left_display_rows, 0)
        self.right_ax.set_ylim(self.right_display_rows, 0)
        self.left_ax.set_xlim(0, self.left_display_cols)
        self.right_ax.set_xlim(0, self.right_display_cols)
        self.left_canvas.draw()
        self.right_canvas.draw()

# ----------------------------------------------------------------
# histograms and stretching
# ----------------------------------------------------------------
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
            left_red_index = self.left_wvl.index(left_red_band)
            left_green_index = self.left_wvl.index(left_green_band)
            left_blue_index = self.left_wvl.index(left_blue_band)

            left_red_channel = np.where(self.left_data[:, :, left_red_index] > 1, np.nan, self.left_data[:, :, left_red_index])
            left_green_channel = np.where(self.left_data[:, :, left_green_index] > 1, np.nan, self.left_data[:, :, left_green_index])
            left_blue_channel = np.where(self.left_data[:, :, left_blue_index] > 1, np.nan, self.left_data[:, :, left_blue_index])

            # Calculate min and max values across all three channels
            left_min_value = min(np.nanmin(left_red_channel), np.nanmin(left_green_channel), np.nanmin(left_blue_channel))
            left_max_value = max(np.nanmax(left_red_channel), np.nanmax(left_green_channel), np.nanmax(left_blue_channel))

            if self.left_hist_window is None or not self.left_hist_window.winfo_exists():
                hist_args = left_min_value, left_max_value, left_red_band, left_green_band, left_blue_band
                self.create_left_histogram(hist_args)

            if self.reset_left_hist:
                red_range = (left_min_value, left_max_value)
                green_range = (left_min_value, left_max_value)
                blue_range = (left_min_value, left_max_value)
                self.reset_left_hist = False
            else:
                red_range =   (self.left_red_min_stretch_var.get() -   abs(self.left_red_min_stretch_var.get()*0.2),   self.left_red_max_stretch_var.get() +   abs(self.left_red_max_stretch_var.get()*0.2))
                green_range = (self.left_green_min_stretch_var.get() - abs(self.left_green_min_stretch_var.get()*0.2), self.left_green_max_stretch_var.get() + abs(self.left_green_max_stretch_var.get()*0.2))
                blue_range =  (self.left_blue_min_stretch_var.get() -  abs(self.left_blue_min_stretch_var.get()*0.2),  self.left_blue_max_stretch_var.get() +  abs(self.left_blue_max_stretch_var.get()*0.2))
            

            # Update the data in the existing subplots
            self.left_hist_axes[0].cla()
            self.left_hist_axes[0].hist(left_red_channel.ravel(), bins=256, color='red', alpha=0.7, range=red_range)
            self.left_hist_axes[0].axvline(self.left_red_min_stretch_var.get(), color='red', linestyle='--', label='Min Stretch')
            self.left_hist_axes[0].axvline(self.left_red_max_stretch_var.get(), color='red', linestyle='--', label='Max Stretch')
            self.left_hist_axes[0].set_title(f'Left Red Channel Histogram: {left_red_band}')
            # self.left_hist_axes[0].legend()

            self.left_hist_axes[1].cla()
            self.left_hist_axes[1].hist(left_green_channel.ravel(), bins=256, color='green', alpha=0.7, range=green_range)
            self.left_hist_axes[1].axvline(self.left_green_min_stretch_var.get(), color='green', linestyle='--', label='Min Stretch')
            self.left_hist_axes[1].axvline(self.left_green_max_stretch_var.get(), color='green', linestyle='--', label='Max Stretch')
            self.left_hist_axes[1].set_title(f'Left Green Channel Histogram: {left_green_band}')
            # self.left_hist_axes[1].legend()

            self.left_hist_axes[2].cla()
            self.left_hist_axes[2].hist(left_blue_channel.ravel(), bins=256, color='blue', alpha=0.7, range=blue_range)
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
            right_red_channel = np.where(self.right_data[:, :, right_red_index] > 1, np.nan, self.right_data[:, :, right_red_index])
            right_green_channel = np.where(self.right_data[:, :, right_green_index] > 1, np.nan, self.right_data[:, :, right_green_index])
            right_blue_channel = np.where(self.right_data[:, :, right_blue_index] > 1, np.nan, self.right_data[:, :, right_blue_index])

            # Calculate min and max values across all three channels
            # right_min_value = min(np.nanmin(right_red_channel), np.nanmin(right_green_channel), np.nanmin(right_blue_channel))
            # right_max_value = max(np.nanmax(right_red_channel), np.nanmax(right_green_channel), np.nanmax(right_blue_channel))

            # calculate the mean and standard deviation of all three channels
            right_red_mean = np.nanmean(right_red_channel)
            right_green_mean = np.nanmean(right_green_channel)
            right_blue_mean = np.nanmean(right_blue_channel)
            right_red_std = np.nanstd(right_red_channel)
            right_green_std = np.nanstd(right_green_channel)
            right_blue_std = np.nanstd(right_blue_channel)

            # set the right_min_value and right_max_value to be -3*std and 3*std from the mean
            right_min_value = min(right_red_mean - 2.5*right_red_std, right_green_mean - 2.5*right_green_std, right_blue_mean - 2.5*right_blue_std)
            right_max_value = max(right_red_mean + 2.5*right_red_std, right_green_mean + 2.5*right_green_std, right_blue_mean + 2.5*right_blue_std)

            if self.right_hist_window is None or not self.right_hist_window.winfo_exists():
                hist_args = right_min_value, right_max_value, right_red_band, right_green_band, right_blue_band
                self.create_right_histogram(hist_args)

            if self.reset_right_hist:
                red_range = (right_min_value, right_max_value)
                green_range = (right_min_value, right_max_value)
                blue_range = (right_min_value, right_max_value)
                self.reset_right_hist = False
            else:
                red_range = (self.right_red_min_stretch_var.get() - abs(self.right_red_min_stretch_var.get()*0.2), self.right_red_max_stretch_var.get() + abs(self.right_red_max_stretch_var.get()*0.2))
                green_range = (self.right_green_min_stretch_var.get() - abs(self.right_green_min_stretch_var.get()*0.2), self.right_green_max_stretch_var.get() + abs(self.right_green_max_stretch_var.get()*0.2))
                blue_range = (self.right_blue_min_stretch_var.get() - abs(self.right_blue_min_stretch_var.get()*0.2), self.right_blue_max_stretch_var.get() + abs(self.right_blue_max_stretch_var.get()*0.2))
            
            # Update the data in the existing subplots
            self.right_hist_axes[0].cla()
            self.right_hist_axes[0].hist(right_red_channel.ravel(), bins=256, color='red', alpha=0.7, range=red_range)
            self.right_hist_axes[0].axvline(self.right_red_min_stretch_var.get(), color='red', linestyle='--', label='Min Stretch')
            self.right_hist_axes[0].axvline(self.right_red_max_stretch_var.get(), color='red', linestyle='--', label='Max Stretch')
            self.right_hist_axes[0].set_title(f'Right Red Channel Histogram: {right_red_band}')
            # self.right_hist_axes[0].legend()

            self.right_hist_axes[1].cla()
            self.right_hist_axes[1].hist(right_green_channel.ravel(), bins=256, color='green', alpha=0.7, range=green_range)
            self.right_hist_axes[1].axvline(self.right_green_min_stretch_var.get(), color='green', linestyle='--', label='Min Stretch')
            self.right_hist_axes[1].axvline(self.right_green_max_stretch_var.get(), color='green', linestyle='--', label='Max Stretch')
            self.right_hist_axes[1].set_title(f'Right Green Channel Histogram: {right_green_band}')
            # self.right_hist_axes[1].legend()

            self.right_hist_axes[2].cla()
            self.right_hist_axes[2].hist(right_blue_channel.ravel(), bins=256, color='blue', alpha=0.7, range=blue_range)
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

    def reset_left_hist_button(self):
        self.reset_left_hist = True
        self.plot_left_histograms()

    def create_left_histogram(self, hist_args):
        left_min_value, left_max_value, left_red_band, left_green_band, left_blue_band = hist_args
        self.left_hist_window = tk.Toplevel(self.root)
        self.left_hist_window.title("Hyperspectral Image Histogram")
        # self.left_hist_window.geometry("600x800")
        
        # Create a frame to hold UI elements with a fixed size
        left_hist_ui_frame = tk.Frame(self.left_hist_window)
        left_hist_ui_frame.pack(fill=tk.X)

        # Create a new histogram plot window if it doesn't exist
        self.left_hist_figure, self.left_hist_axes = plt.subplots(3, 1, figsize=(5, 6))

        self.left_hist_canvas = FigureCanvasTkAgg(self.left_hist_figure, master=self.left_hist_window)
        self.left_hist_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset the histogram
        self.reset_left_hist_x_axis_button = tk.Button(left_hist_ui_frame, text="Reset Histogram", command=self.reset_left_hist_button)
        self.reset_left_hist_x_axis_button.pack(side=tk.RIGHT)

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

    def reset_right_hist_button(self):
        self.reset_right_hist = True
        self.plot_right_histograms()

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
        self.reset_right_hist_x_axis_button = tk.Button(right_hist_ui_frame, text="Reset Histogram", command=self.reset_right_hist_button)
        self.reset_right_hist_x_axis_button.pack(side=tk.RIGHT)

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

# ----------------------------------------------------------------
# polygon functions
# ----------------------------------------------------------------
    def create_polygons_menu_window(self):
        self.polygons_menu_window = tk.Toplevel(self.root)
        self.polygons_menu_window.title("Polygons Menu")
        
        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.polygons_menu_window)
        ui_frame.grid(row=0,column=0, columnspan=7)

        # create variables to track the row and column position of the buttons
        self.polygons_menu_row = 0
        self.polygons_menu_col = 0

        # Create a button to toggle drawing polygons
        self.draw_polygons_button = tk.Button(ui_frame, text="Drawing Polygons Off", command=self.toggle_polygons)
        self.draw_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # create a button to remove polygons from display
        self.remove_polygons_button = tk.Button(ui_frame, text="Remove Polygons from Display", command=self.remove_polygons_from_display)
        self.remove_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # create a button to re-plot the polygons on the display
        self.redraw_all_polygons_button = tk.Button(ui_frame, text="Re-draw Polygons on Display", command=self.draw_all_polygons)
        self.redraw_all_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # create a button to clear all polygons
        self.clear_polygons_button = tk.Button(ui_frame, text="Delete All Polygons", command=self.clear_all_polygons)
        self.clear_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # create a button to extract spectra from polygons
        self.extract_spectra_button = tk.Button(ui_frame, text="Plot Mean Spectra", command=self.update_polygons_spectral_plot)
        self.extract_spectra_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # create a button to save the ROIs
        self.save_rois_button = tk.Button(ui_frame, text="Save ROIs", command=self.save_polygons)
        self.save_rois_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # create a dropdown menu to select the polygon color from a list of colors
        self.polygon_color_var = tk.StringVar(ui_frame)
        self.polygon_color_var.set("red") # default value
        self.polygon_color_dropdown = tk.OptionMenu(ui_frame, self.polygon_color_var, "red", "green", "blue", "yellow", "orange", "purple", "pink", "brown", "black", "white")
        self.polygon_color_dropdown.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_row += 1

        # create a table to display the polygon information
        # # Create headers for the table
        header_labels = ["Polygon Number", "Color", "Number of Points"]
        for i, header_text in enumerate(header_labels):
            header_label = tk.Label(self.polygons_menu_window, text=header_text)
            header_label.grid(row=1, column=i)
        
        # get the information from the polygons to display in the table
        self.update_polygons_table()

    def update_polygons_table(self):
        # get the information from the polygons to display in the table
        self.polygons_table_data = []
        for i, polygon in enumerate(self.all_polygons):
            polygon_number = i
            polygon_color = self.polygon_colors[i]
            number_of_points = len(polygon)-1
            self.polygons_table_data.append([polygon_number, polygon_color, number_of_points])
        
        # Populate the table with data
        for i, row_data in enumerate(self.polygons_table_data):
            for j, cell_value in enumerate(row_data):
                cell_label = tk.Label(self.polygons_menu_window, text=cell_value)
                cell_label.grid(row=i + 2, column=j)

    def clear_polygons_table(self):
        self.polygons_table_data = []
        for i, polygon in enumerate(self.all_polygons):
            polygon_number = None
            polygon_color = None
            number_of_points = None
            self.polygons_table_data.append([polygon_number, polygon_color, number_of_points])
        
        # Populate the table with data
        for i, row_data in enumerate(self.polygons_table_data):
            for j, cell_value in enumerate(row_data):
                cell_label = tk.Label(self.polygons_menu_window, text=cell_value)
                cell_label.grid(row=i + 2, column=j)

    def toggle_polygons(self):
        if self.draw_polygons:
            self.draw_polygons = False
            self.draw_polygons_button.config(text = "Drawing Polygons Off", relief="sunken")
        else:
            self.draw_polygons = True
            self.draw_polygons_button.config(text = "Drawing Polygons On", relief="raised")

    def draw_all_polygons(self):
        for polygon_color, polygon in zip(self.polygon_colors, self.all_polygons):
            x, y = zip(*polygon)
            # self.left_ax.plot(x, y, 'ro')
            # self.right_ax.plot(x, y, 'ro')
            self.left_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color,  edgecolor='k'))
            self.right_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color, edgecolor='k'))
            self.left_canvas.draw()
            self.right_canvas.draw()
            self.update_polygons_table()

    def remove_polygons_from_display(self):
        self.left_ax.clear()
        self.right_ax.clear()
        self.update_left_display()
        self.update_right_display()
    
    def clear_all_polygons(self):
        # Ask the user if they are sure
        answer = messagebox.askyesno("Confirmation", "Are you sure you want to delete all polygons?")
        if answer:
            self.clear_polygons_table()
            self.all_polygons = []
            self.polygon_colors = []
            self.polygon_spectra = []
            self.update_polygons_spectral_plot()
            self.remove_polygons_from_display()
        else:
            pass
    
    def extract_spectra_from_polygons(self):
        if hasattr(self, "left_data"):
            if self.all_polygons:
                # get the indices of the pixels within the polygons
                # left_row_indices = []
                # left_col_indices = []
                self.polygon_spectra = []
                for polygon in self.all_polygons:
                    mask = np.zeros((self.left_data.nrows, self.left_data.ncols), dtype=np.uint8)
                    polygon_points_int = [(int(x), int(y)) for x, y in polygon]
                    cv2.fillPoly(mask, [np.array(polygon_points_int)], 1) 
                    # roi = self.left_data * mask
                    gstats = spectral.calc_stats(self.left_data, mask=mask, allow_nan=True)
                    mean_spectrum = gstats.mean
                    # set values <=0 or >=1 to np.nan
                    mean_spectrum = np.where(mean_spectrum < 0, np.nan, mean_spectrum)
                    mean_spectrum = np.where(mean_spectrum > 1, np.nan, mean_spectrum)

                    self.polygon_spectra.append(mean_spectrum)

                    # polygon = np.array(polygon)
                    # polygon = np.round(polygon).astype(int)
                    # polygon = np.flip(polygon, axis=1)
                    # polygon = np.flip(polygon, axis=0)
                    # polygon = np.transpose(polygon)
                    # row_indices, col_indices = polygon
                    # left_row_indices = np.array(row_indices)
                    # left_col_indices = np.array(col_indices)

                    # # get the spectra from the left data
                    # print(left_row_indices, left_col_indices)
                    # spectrum = self.left_data[left_row_indices, left_col_indices, :]
                    # spectrum = np.where(spectrum < 0, np.nan, spectrum)
                    # spectrum = np.where(spectrum > 1, np.nan, spectrum)
                    # spectrum = np.transpose(spectrum, (1, 0, 2))
                    # spectrum = np.nanmean(spectrum, axis=1)
                    # spectrum = np.nan_to_num(spectrum)

                # plot the spectra
                self.update_polygons_spectral_plot()
            else:
                messagebox.showwarning("Warning", "No polygons drawn. Draw polygons first.")
        else:
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first into left frame.")

    def save_polygons(self):
        if self.all_polygons:
            # get the filename to save the polygons to
            filename = filedialog.asksaveasfilename(initialdir = "/",title = "Select file")#,filetypes = (("shp files","*.shp")))
            if filename:
                # save the polygons using geopandas and shapely.geometry
                all_geo_poly_coords = []
                for poly_coords in self.all_polygons:
                    geo_poly_coords = [self.left_rio.xy(y, x) for x, y in poly_coords]
                    all_geo_poly_coords.append(geo_poly_coords)
                polygon_geometries = [sgp(poly_coords) for poly_coords in all_geo_poly_coords]

                # create a GeoDataFrame
                polygons_gdf = gpd.GeoDataFrame({'color': self.polygon_colors}, geometry=polygon_geometries, crs=self.left_rio.crs)
                polygons_gdf.to_file(filename)
        else:
            messagebox.showwarning("Warning", "No polygons drawn. Draw polygons first.")        

# ----------------------------------------------------------------
# canvas click functionality
# ----------------------------------------------------------------
    def on_left_canvas_click(self, event):
        if hasattr(self, "left_data"):
            if self.draw_polygons:
                if event.inaxes == self.left_ax or event.inaxes == self.right_ax:
                    if event.button == 1:
                        x, y = event.xdata, event.ydata
                        self.current_polygon.append((x, y))
                        point_color = self.polygon_color_var.get()
                        self.left_ax.plot(x, y, color = point_color, marker='o')
                        self.left_canvas.draw()
                        self.right_ax.plot(x, y, color=point_color, marker='o')
                        self.right_canvas.draw()
                    elif event.button == 3:
                        if self.current_polygon:

                            self.current_polygon.append(self.current_polygon[0]) #close the polygon
                            polygon_color = self.polygon_color_var.get()

                            self.left_ax.clear()
                            self.right_ax.clear()
                            self.update_left_display()
                            self.update_right_display()

                            self.all_polygons.append(self.current_polygon)
                            self.polygon_colors.append(polygon_color)
                            self.draw_all_polygons()

                            self.current_polygon = []
                            self.update_polygons_table()
                            self.extract_spectra_from_polygons()


            elif event.inaxes == self.left_ax or event.inaxes == self.right_ax:
                if event.button == 1: #left mouse button press
                    x, y = int(event.xdata), int(event.ydata)
                    if self.bands_inverted:
                        self.spectrum = self.left_data[y, x, :].flatten()
                        self.spectrum = self.spectrum[::-1]
                    else:
                        self.spectrum = self.left_data[y, x, :].flatten()
                    self.spectrum = np.where(self.spectrum < 0, np.nan, self.spectrum)
                    self.spectrum = np.where(self.spectrum > 1, np.nan, self.spectrum)
                    self.update_spectral_plot()
                elif event.button == 3:  # Right mouse button press
                    if hasattr(self, 'spectrum'):
                        # plot ratio spectra
                        x, y = int(event.xdata), int(event.ydata)
                        self.denom_spectrum = self.left_data[y, x, :].flatten()
                        self.denom_spectrum = np.where(self.denom_spectrum < 0, np.nan, self.denom_spectrum)
                        self.denom_spectrum = np.where(self.denom_spectrum > 1, np.nan, self.denom_spectrum)
                        self.ratio_spectrum = self.spectrum / self.denom_spectrum
                        self.update_ratio_spectral_plot()
        else:
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first into left frame.")

    def on_canvas_motion(self, event):
        if hasattr(self, "left_data") and event.inaxes:
            lon, lat = self.left_rio.xy(event.ydata, event.xdata)
            # lon = (event.xdata - self.left_data.ncols/2)*self.x_dpp + float(self.proj_info[3])
            # lat = ((self.left_data.nrows/2) - event.ydata )*self.y_dpp + float(self.proj_info[2])
            coordinate_text = f"Lat: {lat:.4f}, Lon: {lon:.4f}"
            self.coordinates_label.config(text=coordinate_text)
            self.left_canvas.draw()
                
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
                
                self.left_canvas.draw()
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
                
                self.right_canvas.draw()
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
                
                if new_y_lim[0] > self.left_data.nrows:
                    new_y_lim = [self.left_data.nrows, self.left_data.nrows-self.left_display_rows]
                elif new_y_lim[1] < 0:
                    new_y_lim = [self.left_display_rows, 0]
                        
                if new_x_lim[1] > self.left_data.ncols:
                    new_x_lim = [self.left_display_cols-self.left_display_cols, self.left_display_cols]
                elif new_x_lim[0] < 0:
                    new_x_lim = [0, self.left_display_cols]

                self.left_ax.set_xlim(new_x_lim)
                self.left_ax.set_ylim(new_y_lim)
                self.right_ax.set_xlim(new_x_lim)
                self.right_ax.set_ylim(new_y_lim)
                
                self.right_canvas.draw()
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
                
                if new_y_lim[0] > self.left_data.nrows:
                    new_y_lim = [self.left_data.nrows, self.left_data.nrows-self.left_display_rows]
                elif new_y_lim[1] < 0:
                    new_y_lim = [self.left_display_rows, 0]
                        
                if new_x_lim[1] > self.left_data.ncols:
                    new_x_lim = [self.left_display_cols-self.left_display_cols, self.left_display_cols]
                elif new_x_lim[0] < 0:
                    new_x_lim = [0, self.left_display_cols]

                self.left_ax.set_xlim(new_x_lim)
                self.left_ax.set_ylim(new_y_lim)
                self.right_ax.set_xlim(new_x_lim)
                self.right_ax.set_ylim(new_y_lim)
                
                self.left_canvas.draw()
                self.right_canvas.draw()

    def create_left_scroll_bars(self):
        # Create vertical scroll bar
        self.left_vertical_scrollbar = tk.Scrollbar(self.left_frame, orient="vertical", command=self.on_vertical_scroll)
        self.left_vertical_scrollbar.grid(row=0, column=1, rowspan=2, sticky="ns")

        # Create horizontal scroll bar
        self.left_horizontal_scrollbar = tk.Scrollbar(self.left_frame, orient="horizontal", command=self.on_horizontal_scroll)
        self.left_horizontal_scrollbar.grid(row=2, column=0, sticky="ew")

        # Attach scroll bars to the canvas
        self.left_canvas.get_tk_widget().config(yscrollcommand=self.left_vertical_scrollbar.set)
        self.left_canvas.get_tk_widget().config(xscrollcommand=self.left_horizontal_scrollbar.set)

    def create_right_scroll_bars(self):
        # Create vertical scroll bar
        self.right_vertical_scrollbar = tk.Scrollbar(self.right_frame, orient="vertical", command=self.on_vertical_scroll)
        self.right_vertical_scrollbar.grid(row=0, column=1, rowspan=2, sticky="ns")

        # Create horizontal scroll bar
        self.right_horizontal_scrollbar = tk.Scrollbar(self.right_frame, orient="horizontal", command=self.on_horizontal_scroll)
        self.right_horizontal_scrollbar.grid(row=2, column=0, sticky="ew")

        # Attach scroll bars to the canvas
        self.right_canvas.get_tk_widget().config(yscrollcommand=self.right_vertical_scrollbar.set)
        self.right_canvas.get_tk_widget().config(xscrollcommand=self.right_horizontal_scrollbar.set)

    def on_vertical_scroll(self, *args):
        scroll_y = float(args[1])*self.left_display_rows
        y_lim = self.left_ax.get_ylim()
        
        if y_lim[0] + scroll_y > self.left_data.nrows:
            new_lim = [self.left_data.nrows, self.left_data.nrows-(y_lim[0]-y_lim[1])]
        elif y_lim[1] + scroll_y < 0:
            new_lim = [y_lim[0]-y_lim[1], 0]
        else:
            new_lim = [y_lim[0] + scroll_y, y_lim[1] + scroll_y]

        self.left_ax.set_ylim(new_lim[0], new_lim[1])
        self.right_ax.set_ylim(new_lim[0], new_lim[1])

        self.left_canvas.draw()
        self.right_canvas.draw()

        # self.update_left_display()
        # self.update_right_display()

    def on_horizontal_scroll(self, *args):
        scroll_x = float(args[1])*self.left_display_cols
        x_lim = self.left_ax.get_xlim()
        if x_lim[1] + scroll_x > self.left_data.ncols:
            new_lim = [self.left_data.ncols-(x_lim[1]-x_lim[0]), self.left_data.ncols]
        elif x_lim[0] + scroll_x < 0:
            new_lim = [0, x_lim[1] - x_lim[0]]
        else:
            new_lim = [x_lim[0] + scroll_x, x_lim[1] + scroll_x]
        self.left_ax.set_xlim(new_lim[0], new_lim[1])
        self.right_ax.set_xlim(new_lim[0], new_lim[1])

        self.left_canvas.draw()
        self.right_canvas.draw()

# ----------------------------------------------------------------
# Spectral plotting area
# ----------------------------------------------------------------
    # polygons
    # ----------------------------------------------------------------
    def create_polygons_spectral_plot(self):
        self.polygons_spectral_window = tk.Toplevel(self.root)
        self.polygons_spectral_window.title("ROI Spectral Plot")
        
        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.polygons_spectral_window)
        ui_frame.pack(fill=tk.X)

        polygons_spectral_figure, self.polygons_spectral_ax = plt.subplots(figsize=(5,3))
        for poly_color, s in zip(self.polygon_colors, self.polygon_spectra):
            s = s.flatten()
            self.polygon_spectral_lines.append(self.polygons_spectral_ax.plot(self.left_wvl, s, color=poly_color))
        # self.polygon_spectral_lines, = self.polygons_spectral_ax.plot(self.left_wvl, self.polygon_spectra)
        xmin, xmax = self.polygons_spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        all_min_y = []
        all_max_y = []
        for s in self.polygon_spectra:
            min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
            all_min_y.append(min_y)
            all_max_y.append(max_y)
        min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

        self.polygons_spectral_canvas = FigureCanvasTkAgg(polygons_spectral_figure, master=self.polygons_spectral_window)
        self.polygons_spectral_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset x-axis span
        self.reset_polygons_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_polygons_x_axis_span)
        self.reset_polygons_x_axis_button.pack(side=tk.RIGHT)

        # Add built-in span options to a dropdown menu
        span_options = ["Full Span", "410 - 1000 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
        self.polygons_span_var = tk.StringVar()
        self.polygons_span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(ui_frame, textvariable=self.polygons_span_var, values=span_options, state="readonly")
        span_menu.pack(side=tk.RIGHT)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_polygons_x_axis_span)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(self.polygons_spectral_canvas, self.polygons_spectral_window)
        toolbar.update()
        self.polygons_spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Add span selector for x-axis
        self.create_polygons_x_axis_span_selector(self.left_wvl)
    
    def update_polygons_spectral_plot(self):
        if self.polygons_spectral_window is None or not self.polygons_spectral_window.winfo_exists():
            self.polygon_spectral_lines = []
            self.create_polygons_spectral_plot()
        elif self.polygon_spectra:
            # clear everything from the plot before plotting
            self.polygons_spectral_ax.clear()
            self.polygon_spectral_lines = []

            for c, s in zip(self.polygon_colors, self.polygon_spectra):
                s = s.flatten()
                self.polygon_spectral_lines.append(self.polygons_spectral_ax.plot(self.left_wvl, s, color=c))
            xmin, xmax = self.polygons_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            all_min_y = []
            all_max_y = []
            for s in self.polygon_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.polygons_spectral_canvas.draw()
        else:
            self.polygons_spectral_ax.clear()
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
            min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
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
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900)
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
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
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
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)

            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_spectral_canvas.draw()

        self.x_polygons_span_selector = SpanSelector(
            self.polygons_spectral_ax, on_x_span_select, 'horizontal', useblit=True)
    
    # ----------------------------------------------------------------
    # points
    # ----------------------------------------------------------------
    def create_spectral_plot(self):
        self.spectral_window = tk.Toplevel(self.root)
        self.spectral_window.title("Spectral Plot")
        
        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.spectral_window)
        ui_frame.pack(fill=tk.X)

        spectral_figure, self.spectral_ax = plt.subplots(figsize=(5,3))
        self.spectral_line, = self.spectral_ax.plot(self.left_wvl, self.spectrum)
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

        self.spectral_canvas = FigureCanvasTkAgg(spectral_figure, master=self.spectral_window)
        self.spectral_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset x-axis span
        self.reset_spectral_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_x_axis_span)
        self.reset_spectral_x_axis_button.pack(side=tk.RIGHT)

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
        self.create_x_axis_span_selector(self.left_wvl)

    def update_spectral_plot(self):
        if self.spectral_window is None or not self.spectral_window.winfo_exists():
            self.create_spectral_plot()
        else:
            self.spectral_line.set_ydata(self.spectrum)
            xmin, xmax = self.spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.spectral_canvas.draw()

    def reset_x_axis_span(self):
        # Reset x-axis span to the default range
        self.spectral_ax.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
        self.spectral_canvas.draw()

    def update_x_axis_span(self, event):
        selected_span = self.span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_wvl[0], self.left_wvl[-1]),
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
            self.spectral_ax, on_x_span_select, 'horizontal', useblit=True)
    
    # ----------------------------------------------------------------
    # ratio plots
    # ----------------------------------------------------------------
    def create_ratio_spectral_plot(self):
        self.ratio_spectral_window = tk.Toplevel(self.root)
        self.ratio_spectral_window.title("Ratio Spectral Plot")
        
        # Create a frame to hold UI elements with a fixed size
        ratio_ui_frame = tk.Frame(self.ratio_spectral_window)
        ratio_ui_frame.pack(fill=tk.X)

        ratio_spectral_figure, self.ratio_spectral_ax = plt.subplots(figsize=(5,3))
        self.ratio_spectral_line, = self.ratio_spectral_ax.plot(self.left_wvl, self.ratio_spectrum)
        xmin, xmax = self.ratio_spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

        self.ratio_spectral_canvas = FigureCanvasTkAgg(ratio_spectral_figure, master=self.ratio_spectral_window)
        self.ratio_spectral_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Add a button to reset x-axis span
        self.reset_ratio_x_axis_button = tk.Button(ratio_ui_frame, text="Reset X-Axis Span", command=self.reset_ratio_x_axis_span)
        self.reset_ratio_x_axis_button.pack(side=tk.RIGHT)

        # Add built-in span options to a dropdown menu
        span_options = ["Full Span", "410 - 1000 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
        self.ratio_span_var = tk.StringVar()
        self.ratio_span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(ratio_ui_frame, textvariable=self.ratio_span_var, values=span_options, state="readonly")
        span_menu.pack(side=tk.RIGHT)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_ratio_x_axis_span)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(self.ratio_spectral_canvas, self.ratio_spectral_window)
        toolbar.update()
        self.ratio_spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Add span selector for x-axis
        self.create_ratio_x_axis_span_selector(self.left_wvl)

    def update_ratio_spectral_plot(self):
        if self.ratio_spectral_window is None or not self.ratio_spectral_window.winfo_exists():
            self.create_ratio_spectral_plot()
        else:
            self.ratio_spectral_line.set_ydata(self.ratio_spectrum)
            xmin, xmax = self.ratio_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.ratio_spectral_canvas.draw()

    def reset_ratio_x_axis_span(self):
        # Reset x-axis span to the default range
        self.ratio_spectral_ax.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.ratio_spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)
        self.ratio_spectral_canvas.draw()

    def update_ratio_x_axis_span(self, event):
        selected_span = self.ratio_span_var.get()
        # Map the selected span to its corresponding x-axis limits
        span_ranges = {
            "Full Span": (self.left_wvl[0], self.left_wvl[-1]),
            "410 - 1000 nm": (410, 1000),
            "1000 - 2600 nm": (1000, 2600),
            "1200 - 2000 nm": (1200, 2000),
            "1800 - 2500 nm": (1800, 2500),
            "2000 - 2500 nm": (2000, 2500),
            "2700 - 3900 nm": (2700, 3900)
        }
        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            self.ratio_spectral_ax.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))

            # Calculate y-limits based on the data within the new span
            y_min = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx])
            y_max = np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax.set_ylim(y_min - buffer, y_max + buffer)

            self.ratio_spectral_canvas.draw()

    def create_ratio_x_axis_span_selector(self, x_data):
        def on_x_span_select(xmin, xmax):
            xmin_idx = np.argmin(np.abs(x_data - xmin))
            xmax_idx = np.argmin(np.abs(x_data - xmax))
            self.ratio_spectral_ax.set_xlim(x_data[xmin_idx], x_data[xmax_idx])

            # Calculate y-limits based on the data within the new span
            y_min = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx])
            y_max = np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax.set_ylim(y_min - buffer, y_max + buffer)

            self.ratio_spectral_canvas.draw()

        self.ratio_x_span_selector = SpanSelector(
            self.ratio_spectral_ax, on_x_span_select, 'horizontal', useblit=True)

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
                "spectrum": self.spectrum,
                "ratio_spectrum": self.ratio_spectrum,

            }

            # Serialize and save the state dictionary to the chosen file
            with open(file_path, 'wb') as file:
                pickle.dump(state_dict, file)

    def on_closing(self):
        root.destroy()

    def load_state(self):
        # Ask the user to choose the saved session file
        file_path = filedialog.askopenfilename(
            filetypes=[("Pickle files", "*.pkl")],  # Filter for pickle files
        )

        if file_path:
            # Load the saved state from the chosen file
            with open(file_path, 'rb') as file:
                restored_instance = pickle.load(file)

            # Update the current instance with the restored state
            self.load_left_data(restored_instance['left_file'])
            self.load_right_data(restored_instance['right_file'])
            self.all_polygons = restored_instance['all_polygons']
            self.polygon_colors = restored_instance['polygon_colors']
            self.polygon_spectra = restored_instance['polygon_spectra']
            self.spectrum = restored_instance['spectrum']
            self.ratio_spectrum = restored_instance['ratio_spectrum']

if __name__ == "__main__":
    root = tk.Tk()
    app = SpectralCubeAnalysisTool(root)
    root.protocol("WM_DELETE_WINDOW", app.on_closing)
    root.mainloop()
