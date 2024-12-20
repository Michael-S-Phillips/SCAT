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
import cv2
import geopandas as gpd
from shapely.geometry import Polygon as sgp
from scipy import interpolate
import pandas as pd
import glob
import matplotlib.colors as mcolors
from joblib import Parallel, delayed
import os

# try:
from hypyrameter.paramCalculator import cubeParamCalculator
from hypyrameter.interpNans import interpNaNs as interpNaNs
from hypyrameter.utils import getSpecFiles, getWavelengthFromUSGS, getReflectanceFromUSGS
# except:
#     print("Unable to load HyPyRameter, visit https://github.com/Michael-S-Phillips/HyPyRameter for more information.")

class SpectralCubeAnalysisTool:
    '''
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
    '''
# ----------------------------------------------------------------
# initial setup
# ----------------------------------------------------------------
    def __init__(self, root):
        '''
        Initializes the GUI and creates the main UI elements.
        '''
        self.root = root
        self.root.title("Spectral Cube Analysis Tool")

        self.default_rgb_bands = [3, 2, 1]  # Default bands for RGB display of spectral data
        self.default_parameter_bands = [3, 2, 1]  # Default bands for RGB display of parameter data, MAF
        self.default_stretch = 'linear'  # Default stretch for RGB display
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
        self.left_is_paramenter = False
        self.reset_right_hist = False
        self.reset_left_hist = False
        self.draw_polygons = False
        self.ignore_bad_bands_flag = False
        self.ignore_polygons_bad_bands_flag = False
        self.ignore_ratio_bad_bands_flag = False
        self.draw_poly_motion_flag = False
        self.lock_column = False
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
        self.collected_points = pd.DataFrame(columns=("name", "color", "mineral_id1", "mineral_id2", "mineral_id3", "mineral_id4", "x", "y", "lon", "lat", "spectrum", "denom", "parameter_values"))
        self.polygons_on_flag = False
        self.upload_polygons_flag = False

        self.usgs_spectra_path = 'librarySpectra/'
        self.usgs_spectra_folders = glob.glob(self.usgs_spectra_path+'/*')
        self.mica_spectra_path = 'mrocr_8001/mrocr_8001/data/'
        self.mica_spectra_names = glob.glob(self.mica_spectra_path+'*.tab')

        self.pc = None

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

        # Add a button to calculate spectral parameters
        # self.parameter_calculation_button = tk.Button(self.right_buttons_frame, text="Calculate Spectral Parameters", command=self.calculate_spectral_parameters)
        # self.parameter_calculation_button.grid(row=2, column=0, sticky = 'new')

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

        self.apply_button = tk.Button(self.button_frame, text="Apply", command=self.apply_new_left_bands)
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

        self.apply_button = tk.Button(self.button_frame, text="Apply", command=self.apply_new_right_bands)
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
        file_menu.add_command(label = "Save Session", command=self.save_state)
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        histogram_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Plot Hitograms", menu=histogram_menu)
        histogram_menu.add_command(label="Plot Left Frame Histogram", command=self.plot_left_histograms)
        histogram_menu.add_command(label="Plot Right Frame Histogram", command=self.plot_right_histograms)

        processing_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Processing", menu=processing_menu)
        processing_menu.add_command(label="Calculate Spectral Parameters", command=self.calculate_spectral_parameters)
        processing_menu.add_command(label="Select Bad Bands", command=self.select_bad_bands)

    def create_left_canvas(self):
        self.left_figure, self.left_ax = plt.subplots(figsize=(6, 6))
        self.left_ax.set_title(f"Hyperspectral Cube: {self.img_base_name}", fontsize= 8)
        self.left_ax.tick_params(axis='both', which='major', labelsize=8)
        # Set the margins and spacing to zero
        self.left_figure.tight_layout()
        self.left_figure.subplots_adjust(wspace=0, hspace=0)
        
        self.left_frame = tk.Frame(self.main_ui_frame)
        self.left_frame.grid(row=1, column=0, sticky="nsew")

        self.left_canvas = FigureCanvasTkAgg(self.left_figure, master=self.left_frame)
        self.left_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")  # Use grid instead of pack

        self.left_nav_toolbar_frame = tk.Frame(self.left_frame)
        self.left_nav_toolbar = NavigationToolbar2Tk(self.left_canvas, self.left_nav_toolbar_frame)
        self.left_nav_toolbar.update()
        self.left_nav_toolbar_frame.grid(row=1, column=0, sticky="nsew")

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
        self.left_canvas.mpl_connect('button_release_event', self.on_left_release)
        # self.left_canvas.mpl_connect('draw_event', self.on_scroll)
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

        self.right_nav_toolbar_frame = tk.Frame(self.right_frame)
        self.right_nav_toolbar = NavigationToolbar2Tk(self.right_canvas, self.right_nav_toolbar_frame)
        self.right_nav_toolbar.update()
        self.right_nav_toolbar_frame.grid(row=1, column=0, sticky="nsew")


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
        self.right_canvas.mpl_connect('button_release_event', self.on_right_release)
        # self.right_canvas.mpl_connect('button_release_event', self.on_canvas_release)
        self.right_canvas.mpl_connect('motion_notify_event', self.on_canvas_motion)

# ----------------------------------------------------------------
# loading the data
# ----------------------------------------------------------------
    def load_left_data(self, file_path = None):
        if not file_path:
            file_path = filedialog.askopenfilename(filetypes=[("ENVI Files", "*.hdr"), ("JP2 Files", "*.jp2"), ("JP2 Files (Uppercase)", "*.JP2")])
            # file_ext = file_path.split('.')[-1]
            # print(file_ext, file_path)
        if file_path:
            # if file_ext == 'hdr' or file_ext == 'HDR':
            self.left_data = spectral.io.envi.open(file_path)
            self.ld = self.left_data.load()
            self.img_base_name = self.left_data.filename.split('/')[-1].split('.')[0]
            rio_path = file_path.replace("hdr", "img")
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

            name_ = file_path.split('/')[-1].split('.')[0]
            self.root.title(f"Spectral Cube Analysis Tool: {name_}")

    def load_right_data(self, file_path = None):
        if not file_path:
            file_path = filedialog.askopenfilename(filetypes=[("ENVI Files", "*.hdr"), ("JP2 Files", "*.jp2"), ("JP2 Files (Uppercase)", "*.JP2")])
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
                        if self.left_data.bands.centers[1] - self.left_data.bands.centers[-1] > 0:
                            self.bands_inverted = True
                        if self.bands_inverted:
                            wavelengths = self.left_data.bands.centers
                            self.left_wvl = sorted([float(i) for i in self.left_data.bands.centers])
                        else:
                            wavelengths = self.left_data.bands.centers
                            self.left_wvl = [float(i) for i in self.left_data.bands.centers]

                elif self.left_data.metadata['band names'] is not None:
                    if str(self.left_data.metadata['band names'][1]).isalpha():
                        # This is a band parameter image
                        self.left_is_parameter = True
                        wavelengths = self.left_data.metadata['band names']
                        self.left_wvl = self.left_data.metadata['band names']
                    else:
                        # This is a spectral cube
                        wavelengths = np.array(self.left_data.metadata['band names']).astype(float).tolist()
                        self.left_data.metadata['band names'] = wavelengths
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
                wavelengths = list(range(1, self.left_data.shape[2] + 1))
                self.left_wvl = wavelengths
                messagebox.showerror("Error", "Unable to load wavelength information. Setting wavelengths to default integer values")
            import numpy as np

            # set the wavelength values
            self.left_red_band_menu["values"] = wavelengths
            self.left_green_band_menu["values"] = wavelengths
            self.left_blue_band_menu["values"] = wavelengths
            self.left_red_band_menu.set(wavelengths[self.default_rgb_bands[0]])
            self.left_green_band_menu.set(wavelengths[self.default_rgb_bands[1]])
            self.left_blue_band_menu.set(wavelengths[self.default_rgb_bands[2]])

            def calculate_stretch_limits(image, channel_index):
                # Calculate outlier before defining stretch limits
                channel_data = image[:, :, channel_index]
                channel_data = np.where((channel_data > -1) & (channel_data < 1), channel_data, np.nan)
                mean_value = np.nanmean(channel_data)
                std_dev = np.nanstd(channel_data)
                z_scores = (channel_data - mean_value) / std_dev
                
                # Define a threshold for considering values as outliers (e.g., z-score > 3)
                outliers = np.abs(z_scores) > 3
                
                # Calculate statistics on non-outlier values
                non_outliers = channel_data[~outliers]
                
                # Calculate stretch limits
                min_stretch = np.nanmean(non_outliers) - 2*np.nanstd(non_outliers)
                max_stretch = np.nanmean(non_outliers) + 2*np.nanstd(non_outliers)
                
                return min_stretch, max_stretch

            # default_rgb_bands contains the indices for red, green, and blue channels
            red_index, green_index, blue_index = self.default_rgb_bands

            # Calculate stretch limits for each channel
            red_min_stretch, red_max_stretch = calculate_stretch_limits(self.left_data, red_index)
            green_min_stretch, green_max_stretch = calculate_stretch_limits(self.left_data, green_index)
            blue_min_stretch, blue_max_stretch = calculate_stretch_limits(self.left_data, blue_index)

            # Set the stretch limits variables for each channel
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
                wavelengths = list(range(1, self.right_data.shape[2] + 1))
                self.right_wvl = wavelengths
                messagebox.showerror("Error", "Unable to load wavelength information.")

            if len(wavelengths) < np.max(self.default_parameter_bands):
                self.default_parameter_bands = [2, 1, 0]

            self.right_red_band_menu["values"] = wavelengths
            self.right_green_band_menu["values"] = wavelengths
            self.right_blue_band_menu["values"] = wavelengths
            self.right_red_band_menu.set(wavelengths[self.default_parameter_bands[0]])
            self.right_green_band_menu.set(wavelengths[self.default_parameter_bands[1]])
            self.right_blue_band_menu.set(wavelengths[self.default_parameter_bands[2]])

            def calculate_stretch_limits(image, channel_index):
                # Calculate median and median absolute deviation (MAD)
                channel_data = image[:, :, channel_index]
                median_value = np.nanmedian(channel_data)
                mad = np.nanmedian(np.abs(channel_data - median_value))

                # Define a threshold for considering values as outliers (e.g., 3 times MAD)
                threshold = 2 * mad
                
                # Clip values beyond the threshold to focus on the main Gaussian distribution
                channel_data_clipped = np.clip(channel_data, median_value - threshold, median_value + threshold)

                # Calculate stretch limits
                min_stretch = np.nanmean(channel_data_clipped) - 2*np.nanstd(channel_data_clipped)
                max_stretch = np.nanmean(channel_data_clipped) + 2*np.nanstd(channel_data_clipped)

                return min_stretch, max_stretch

            # default_rgb_bands contains the indices for red, green, and blue channels
            red_index, green_index, blue_index = self.default_parameter_bands

            # Calculate stretch limits for each channel in the "right_data"
            right_red_min_stretch, right_red_max_stretch = calculate_stretch_limits(self.right_data, red_index)
            right_green_min_stretch, right_green_max_stretch = calculate_stretch_limits(self.right_data, green_index)
            right_blue_min_stretch, right_blue_max_stretch = calculate_stretch_limits(self.right_data, blue_index)

            # Set the stretch limits variables for each channel in the "right_data"
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

            if self.left_hist_window is not None and self.left_hist_window.winfo_exists():
                self.plot_left_histograms()

        except ValueError:
            messagebox.showerror(f"Error", "Invalid wavelength. Please enter valid wavelengths.\n{ValueError}")

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

            if self.right_hist_window is not None and self.right_hist_window.winfo_exists():
                self.plot_right_histograms()

        except ValueError:
            messagebox.showerror(f"Error", "Invalid wavelength. Please enter valid wavelengths.\n{ValueError}")

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
        self.left_canvas.draw()
        # show(np.transpose(left_rgb_image, (2,0,1)), ax=self.left_ax, transform=self.left_transform)

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

    def apply_new_left_bands(self):
        red_band = float(self.left_red_band_var.get())
        green_band = float(self.left_green_band_var.get())
        blue_band = float(self.left_blue_band_var.get())
        red_index = self.left_wvl.index(red_band)
        green_index = self.left_wvl.index(green_band)
        blue_index = self.left_wvl.index(blue_band)

        self.left_red_min_stretch_var.set(np.nanmedian(self.left_data[:,:,red_index])-np.nanquantile(self.left_data[:,:,red_index], 0.6))
        self.left_red_max_stretch_var.set(np.nanmedian(self.left_data[:,:,red_index])+np.nanquantile(self.left_data[:,:,red_index], 0.6))
        self.left_green_min_stretch_var.set(np.nanmedian(self.left_data[:,:,green_index])-np.nanquantile(self.left_data[:,:,green_index], 0.6))
        self.left_green_max_stretch_var.set(np.nanmedian(self.left_data[:,:,green_index])+np.nanquantile(self.left_data[:,:,green_index], 0.6))
        self.left_blue_min_stretch_var.set(np.nanmedian(self.left_data[:,:,blue_index])-np.nanquantile(self.left_data[:,:,blue_index], 0.6))
        self.left_blue_max_stretch_var.set(np.nanmedian(self.left_data[:,:,blue_index])+np.nanquantile(self.left_data[:,:,blue_index], 0.6))

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

        self.right_red_min_stretch_var.set(np.nanmedian(self.right_data[:,:,red_index])-np.nanquantile(self.right_data[:,:,red_index], 0.6))
        self.right_red_max_stretch_var.set(np.nanmedian(self.right_data[:,:,red_index])+np.nanquantile(self.right_data[:,:,red_index], 0.6))
        self.right_green_min_stretch_var.set(np.nanmedian(self.right_data[:,:,green_index])-np.nanquantile(self.right_data[:,:,green_index], 0.6))
        self.right_green_max_stretch_var.set(np.nanmedian(self.right_data[:,:,green_index])+np.nanquantile(self.right_data[:,:,green_index], 0.6))
        self.right_blue_min_stretch_var.set(np.nanmedian(self.right_data[:,:,blue_index])-np.nanquantile(self.right_data[:,:,blue_index], 0.6))
        self.right_blue_max_stretch_var.set(np.nanmedian(self.right_data[:,:,blue_index])+np.nanquantile(self.right_data[:,:,blue_index], 0.6))

        self.update_right_display()

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
            right_min_value = min(right_red_mean - right_red_std, right_green_mean - right_green_std, right_blue_mean - right_blue_std)
            right_max_value = max(right_red_mean + right_red_std, right_green_mean + right_green_std, right_blue_mean + right_blue_std)

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
            if self.draw_polygons and not self.left_nav_toolbar.mode == 'pan/zoom' and not self.right_nav_toolbar.mode == 'pan/zoom':
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
                            if self.template_index == self.template_index_cache and self.current_denom_polygon_index is not None:
                                print(self.template_index, self.template_index_cache)
                                # return the item from self.polygon_table at row index self.current_denom_polygon_index
                                # if self.current_denom_polygon_index is None:
                                #     self.current_denom_polygon_index = len(self.all_polygons)
                                item = self.polygon_table_item_id_row[self.current_denom_polygon_index]
                                self.polygon_table.delete(item)
                                self.all_polygons.pop(self.current_denom_polygon_index)
                                self.polygon_colors.pop(self.current_denom_polygon_index)
                                self.best_denom_id_list.pop(self.current_denom_polygon_index)
                                self.ratio_spectrum_list.pop(self.current_denom_polygon_index)
                                if self.polygon_spectra:
                                    self.polygon_spectra.pop(self.current_denom_polygon_index)
                                if self.polygons_spectral_window is not None:
                                    if self.polygons_spectral_window.winfo_exists():
                                        self.update_polygons_spectral_plot()
                                self.remove_polygons_from_display()
                                self.draw_all_polygons()

                            self.template_index_cache = self.template_index
                            template_polygon = self.all_polygons[self.template_index]
                            denom_x_centroid, denom_y_centroid = event.xdata, event.ydata
                            x_centroid, y_centroid = np.nanmean(template_polygon, axis=0)
                            dx, dy = denom_x_centroid-x_centroid, denom_y_centroid-y_centroid
                            if self.lock_column:
                                dx = 0
                            self.denom_polygon = [(point[0]+dx, point[1]+dy) for point in template_polygon]

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


            elif event.inaxes == self.left_ax or event.inaxes == self.right_ax and not self.left_nav_toolbar.mode == 'pan/zoom' and not self.right_nav_toolbar.mode == 'pan/zoom':
                if event.button == 1: #left mouse button press
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
                    self.left_ax.plot(x, y, 'kx')
                    self.right_ax.plot(x, y, 'kx')
                    self.left_canvas.draw()
                    self.right_canvas.draw()
                elif event.button == 3:  # Right mouse button press
                    if hasattr(self, 'spectrum'):
                        x, y = int(event.xdata), int(event.ydata)
                    
                        # add an x on the image frames where the spectrum was clicked
                        self.left_ax.plot(x, y, 'rx')
                        self.right_ax.plot(x, y, 'rx')
                        self.left_canvas.draw()
                        self.right_canvas.draw()

                        # plot ratio spectra
                        self.denom_spectrum = self.left_data[y, x, :].flatten()
                        self.denom_spectrum = np.where(self.denom_spectrum < 0, np.nan, self.denom_spectrum)
                        self.denom_spectrum = np.where(self.denom_spectrum > 1.5, np.nan, self.denom_spectrum)
                        self.ratio_spectrum = self.spectrum / self.denom_spectrum
                        self.update_ratio_spectral_plot()
        else:
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first into left frame.")

    def on_canvas_motion(self, event):
        if hasattr(self, "left_data") and event.inaxes:
            lon, lat = self.left_rio.xy(event.ydata, event.xdata)
            coordinate_text = f"Lat: {lat:.4f}, Lon: {lon:.4f}"
            self.coordinates_label.config(text=coordinate_text)
            self.left_canvas.draw()
        if self.draw_polygons and self.draw_poly_motion_flag and event.inaxes and not self.left_nav_toolbar.mode == 'pan/zoom':
            x, y = event.xdata, event.ydata
            self.current_polygon.append((x, y))
            point_color = self.polygon_color_var.get()
            self.left_ax.plot(x, y, color = point_color, marker='.')
            self.left_canvas.draw()
            self.right_ax.plot(x, y, color=point_color, marker='.')
            self.right_canvas.draw()
        
    def on_left_release(self, event):
        if hasattr(self, "left_data") and self.left_ax.in_axes(event):
            if self.left_nav_toolbar.mode == 'pan/zoom':
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
                    self.current_polygon.append(self.current_polygon[0]) #close the polygon
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
                    messagebox.showwarning("Warning", "Polygons must have more than 2 points.")
    
    def on_right_release(self, event):
        if hasattr(self, "right_data") and self.right_ax.in_axes(event):
            if self.right_nav_toolbar.mode == 'pan/zoom':
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
                    self.current_polygon.append(self.current_polygon[0]) #close the polygon
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
                    messagebox.showwarning("Warning", "Polygons must have more than 2 points.")
                
    def on_scroll(self, event):
        # Handle scroll event for the left canvas
        x_data, y_data = event.xdata, event.ydata
        x_lim, y_lim = self.left_ax.get_xlim(), self.left_ax.get_ylim()
        
        # Define the zoom factor (adjust as needed)
        zoom_factor = 1.1 if event.button == 'down' else 1 / 1.1  # Zoom in or out

        # Adjust the axis limits centered on the mouse cursor position
        new_x_lim = [x_data - (x_data - x_lim[0]) * zoom_factor, x_data + (x_lim[1] - x_data) * zoom_factor]
        new_y_lim = [y_data - (y_data - y_lim[0]) * zoom_factor, y_data + (y_lim[1] - y_data) * zoom_factor]

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
        ui_frame.grid(row=0,column=0, columnspan=7)

        # ----------------------------------------------------------------
        # create variables to track the row and column position of the buttons
        self.polygons_menu_row = 0
        self.polygons_menu_col = 0

        # ----------------------------------------------------------------
        # Create a button to toggle drawing polygons
        self.draw_polygons_button = tk.Button(ui_frame, text="Drawing Polygons Off", command=self.toggle_polygons)
        self.draw_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to toggle column lock
        self.lock_column_button = tk.Button(ui_frame, text="Column Lock Off", command=self.toggle_column_lock)
        self.lock_column_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to remove polygons from display
        self.remove_polygons_button = tk.Button(ui_frame, text="Remove Polygons from Display", command=self.remove_polygons_from_display)
        self.remove_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to re-plot the polygons on the display
        self.redraw_all_polygons_button = tk.Button(ui_frame, text="Re-draw Polygons on Display", command=self.draw_all_polygons)
        self.redraw_all_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to clear all polygons
        self.clear_polygons_button = tk.Button(ui_frame, text="Delete All Polygons", command=self.clear_all_polygons)
        self.clear_polygons_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to extract spectra from polygons
        self.extract_spectra_button = tk.Button(ui_frame, text="Plot Mean Spectra", command=self.update_polygons_spectral_plot)
        self.extract_spectra_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a button to save the ROIs
        self.save_rois_button = tk.Button(ui_frame, text="Save ROIs", command=self.save_polygons)
        self.save_rois_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1
        
        # ----------------------------------------------------------------
        # create a button to load the ROIs
        self.load_rois_button = tk.Button(ui_frame, text="Load ROIs", command=self.load_polygons)
        self.load_rois_button.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_col += 1

        # ----------------------------------------------------------------
        # create a dropdown menu to select the polygon color from a list of colors
        self.polygon_color_var = tk.StringVar(ui_frame)
        self.polygon_color_var.set("black") # default value
        mcolors_keys = mcolors.CSS4_COLORS.keys()
        colors_list = [value for value in mcolors_keys]
        self.polygon_color_dropdown = tk.OptionMenu(ui_frame, self.polygon_color_var, *colors_list)
        # self.polygon_color_dropdown = tk.OptionMenu(ui_frame, self.polygon_color_var, "red", "green", "blue", "yellow", "orange", "purple", "pink", "brown", "black", "white")
        self.polygon_color_dropdown.grid(row=self.polygons_menu_row, column=self.polygons_menu_col)
        self.polygons_menu_row += 1

        # ----------------------------------------------------------------
        # create a table to display polygon information
        self.create_polygons_table()

    # ------------------------------------------------
    # polygon table functions
    # ------------------------------------------------
    def create_polygons_table(self):
        self.polygon_table_header_labels = ("Polygon Number", "Color", "Number of Points", "Denominator", "Template", "Mineral ID 1", "Mineral ID 2", "Mineral ID 3", "Mineral ID 4", "wvl", "Spectrum Mean", "params", "Parameters Mean", "Best Denom ID", "Ratio Spectrum")
        self.polygon_editing_data = {} 
        self.editing_polygon_table = False

        # Create a Text widget for displaying the header labels
        self.polygon_table = ttk.Treeview(self.polygons_menu_window, columns=self.polygon_table_header_labels, show="headings")

        # Configure column headings
        for col in self.polygon_table_header_labels:
            self.polygon_table.heading(col, text=col)
            self.polygon_table.column(col, width=len(col)*9, anchor=tk.CENTER)  # Set the column width as needed
        
        self.polygon_table.grid(row=1, column=0, columnspan=len(self.polygon_table_header_labels))

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
        self.context_menu.add_command(label="Plot Spectrum", command=self.plot_selected_polygon)
        self.context_menu.add_command(label="Plot Ratio Spectrum", command=lambda: self.plot_selected_polygon(d = True))

        self.update_polygons_table()

    def update_polygons_table(self):
        # get the information from the polygons to display in the table 13 14 
        if hasattr(self, 'polygons_table_data'):
            if len(self.polygons_table_data) != 0:
                min_ids = [line[5:9] for line in self.polygons_table_data]
                if self.upload_polygons_flag == False:
                    try:
                        self.best_denom_id_list = [line[13] for line in self.polygons_table_data]
                        self.ratio_spectrum_list = [line[14] for line in self.polygons_table_data]
                    except:
                        print('No best denom id or ratio spectrum')
        self.polygons_table_data = []
        for i, polygon in enumerate(self.all_polygons):
            polygon_number = i
            polygon_color = self.polygon_colors[i]
            number_of_points = len(polygon) - 1
            denom_check = '[ ]'
            template_check = '[ ]'
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
                denom_check = '[X]'
            if self.template_index == i:
                template_check = '[X]'
            if 'min_ids' in locals() and i<len(min_ids):
                minid1, minid2, minid3, minid4 = min_ids[i]
            else:
                minid1, minid2, minid3, minid4 = '', '', '', ''
            try:
                spec_mean = self.polygon_spectra[i]
                param_mean = self.polygon_params[i]
            except IndexError:
                print(f'No means for polygon {i}')
            
            row = (polygon_number, polygon_color, number_of_points, denom_check, template_check, minid1, minid2, minid3, minid4, self.left_wvl, spec_mean, self.right_wvl, param_mean, best_denom_id, row_ratio_spectrum)
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
            if cell_value == '[ ]':
                # set all the other items in the column to '[ ]'
                for item_ in self.polygon_table.get_children():
                    self.polygon_table.set(item_, column, '[ ]')
                self.polygon_table.set(item, column, '[X]')
                # populate self.denominator index with the value of polygon Number column
                self.denominator_index = int(self.polygon_table.item(item, "values")[0])
            else:
                self.polygon_table.set(item, column, '[ ]')
                self.denominator_index = None

        if int(column[1:]) == 5:
            if cell_value == '[ ]':
                # set all the other items in the column to '[ ]'
                for item_ in self.polygon_table.get_children():
                    self.polygon_table.set(item_, column, '[ ]')
                self.polygon_table.set(item, column, '[X]')
                # populate self.denominator index with the value of polygon Number column
                self.template_index = int(self.polygon_table.item(item, "values")[0])
            else:
                self.polygon_table.set(item, column, '[ ]')
                self.template_index = []
                self.template_index_cache = []
        
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
   
    def show_context_menu(self, event):
        self.tmp_event = event
        item = self.polygon_table.identify_row(event.y)
        if item:
            self.polygon_table.selection_set(item)
            self.context_menu.post(event.x_root, event.y_root)

    def edit_cell(self, event=None):
        # if event is None:
        #     event = self.tmp_event
        if event is not None:
            self.tmp_event = event
        event = self.tmp_event

        item = self.polygon_table.identify_row(event.y)
        column = self.polygon_table.identify_column(event.x)
        if int(column[1:]) not in (2, 6, 7, 8, 9):
            return
        cell_value = self.polygon_table.item(item, "values")[int(column[1:]) - 1]

        if not self.editing_polygon_table:
            self.polygon_editing_data[item, column] = cell_value  # Store the original value

            # Create a Toplevel window for editing
            edit_window = tk.Toplevel(self.polygons_menu_window)

            # Create an Entry widget inside the edit_window
            edit_entry = ttk.Entry(edit_window, justify="center")
            edit_entry.insert(0, cell_value)
            edit_entry.pack()

            def save_and_close():
                new_value = edit_entry.get()
                self.polygon_table.set(item, column, new_value)
                pc_index = int(self.polygon_table.item(item, "values")[0])
                if int(column[1:]) == 2:
                    self.polygon_colors[pc_index] = new_value
                else:
                    updated_row = list(self.polygons_table_data[pc_index])
                    updated_row[int(column[1:]) - 1] = new_value
                    self.polygons_table_data[pc_index] = tuple(updated_row)
                self.draw_all_polygons()
                self.update_polygons_spectral_plot()
                edit_window.destroy()

            # Create a button to save changes
            save_button = ttk.Button(edit_window, text="Save", command=save_and_close)
            save_button.pack()

            # Bind the "Return" key to the "Save" button's functionality
            edit_window.bind("<Return>", lambda event: save_and_close())

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
            self.draw_polygons_button.config(text = "Drawing Polygons Off", relief="sunken")
        else:
            self.draw_polygons = True
            self.draw_polygons_button.config(text = "Drawing Polygons On", relief="raised")

    def toggle_column_lock(self):
        if self.lock_column:
            self.lock_column = False
            self.lock_column_button.config(text = "Lock Column Off", relief="sunken")
        else:
            self.lock_column = True
            self.lock_column_button.config(text = "Lock Column On", relief="raised")
    
    def draw_all_polygons(self):
        self.polygons_on_flag = True
        for polygon_color, polygon in zip(self.polygon_colors, self.all_polygons):
            # x, y = zip(*polygon)
            # self.left_ax.plot(x, y, 'ro')
            # self.right_ax.plot(x, y, 'ro')
            if self.highlight_polygon:
                if self.all_polygons.index(polygon) == self.polygon_to_highlight:
                    self.left_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color,  edgecolor='cyan'))
                    self.right_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color, edgecolor='cyan'))
                else:
                    self.left_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color,  edgecolor='k'))
                    self.right_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color, edgecolor='k'))
            else:
                self.left_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color,  edgecolor='k'))
                self.right_ax.add_patch(Polygon(polygon, closed=True, facecolor=polygon_color, edgecolor='k'))
        self.left_canvas.draw()
        self.right_canvas.draw()
        self.update_polygons_table()

    def remove_polygons_from_display(self):
        self.polygons_on_flag = False
        self.clear_axes()
        self.update_left_display()
        self.update_right_display()
    
    def delete_polygon(self):
        event = self.tmp_event
        # Ask the user if they are sure
        answer = messagebox.askyesno("Confirmation", "Are you sure you want to delete this polygon?")
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

    def plot_selected_polygon(self, d = None):
        event = self.tmp_event
        item = self.polygon_table.identify_row(event.y)
        pc_index = int(self.polygon_table.item(item, "values")[0])
        ssw = spectral_window #single_spectrum_window
        # ssw.create_window(app, pc_index) #old
        ssw.create_spectral_plot(app, pc_index, d = d)
        # spec = self.polygon_spectra[pc_index]

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

            # Function to calculate statistics for a single mask
            def calculate_stats(data, mask, allow_nan):
                gstats = spectral.calc_stats(data, mask=mask, allow_nan=allow_nan)
                mean_spectrum = gstats.mean
                mean_spectrum = np.where(mean_spectrum < 0, np.nan, mean_spectrum)
                mean_spectrum = np.where(mean_spectrum > 1.5, np.nan, mean_spectrum)
                return mean_spectrum

            mean_spec = calculate_stats(self.ld, mask, True)
            mean_param = calculate_stats(self.rd, mask, True)
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
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first into left frame.")

    def extract_spectra_from_polygons(self):
        if hasattr(self, "left_data"):
            self.polygon_spectra = []
            self.polygon_params = []

            if self.all_polygons:
                nrows, ncols = self.left_data.nrows, self.left_data.ncols

                # Prepare masks for all polygons
                print(f'Extracting means for {len(self.all_polygons)} polygons')
                masks = np.zeros((len(self.all_polygons), nrows, ncols), dtype=np.uint8)
                for idx, polygon in enumerate(self.all_polygons):
                    polygon_points_int = [(int(x), int(y)) for x, y in polygon]
                    cv2.fillPoly(masks[idx], [np.array(polygon_points_int)], 1)

                # Function to calculate statistics for a single mask
                def calculate_stats(data, mask, allow_nan):
                    gstats = spectral.calc_stats(data, mask=mask, allow_nan=allow_nan)
                    mean_spectrum = gstats.mean
                    mean_spectrum = np.where(mean_spectrum < 0, np.nan, mean_spectrum)
                    mean_spectrum = np.where(mean_spectrum > 1.5, np.nan, mean_spectrum)
                    return mean_spectrum

                # Calculate statistics in parallel
                results = Parallel(n_jobs=-1)(delayed(calculate_stats)(self.ld, mask, True) for mask in masks)
                print(len(results))
                self.polygon_spectra.extend(results)

                # Calculate parameters (assuming self.right_data is defined)
                results_params = Parallel(n_jobs=-1)(delayed(calculate_stats)(self.rd, mask, True) for mask in masks)
                self.polygon_params.extend(results_params)

                # plot the spectra
                self.update_polygons_spectral_plot()
            else:
                messagebox.showwarning("Warning", "No polygons drawn. Draw polygons first.")
        else:
            messagebox.showwarning("Warning", "No data loaded. Load hyperspectral data first into left frame.")


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
            filename = filedialog.asksaveasfilename(initialdir = "/",title = "Select file")#,filetypes = (("shp files","*.shp")))
            if filename:
                # Add .gpkg extension if not already present
                if not filename.lower().endswith('.gpkg'):
                    filename += '.gpkg'

                # if the file is already present, ask the user if they want to save
                # Check if file already exists
                if os.path.exists(filename):
                    # Ask user if they want to overwrite the file
                    overwrite = messagebox.askyesno("Overwrite Confirmation", f"The file '{filename}' already exists. Do you want to overwrite it?")
                    if not overwrite:
                        return
                # save the polygons using geopandas and shapely.geometry
                all_geo_poly_coords = []
                for poly_coords in self.all_polygons:
                    geo_poly_coords = [self.left_rio.xy(y, x) for x, y in poly_coords]
                    all_geo_poly_coords.append(geo_poly_coords)
                polygon_geometries = [sgp(poly_coords) for poly_coords in all_geo_poly_coords]

                # create a GeoDataFrame     
                column_dict = {}
                # Assuming each column has a unique identifier (replace 'column_id' with the actual identifier)
                for i, column_id in enumerate(self.polygon_table_header_labels):  # Add the actual identifiers for your columns
                    column_values = [self.polygon_table.item(item_id, 'values')[i] for item_id in self.polygon_table.get_children()]
                    column_dict[column_id] = column_values
                print(column_dict)
                polygons_gdf = gpd.GeoDataFrame(column_dict, geometry=polygon_geometries, crs=self.left_rio.crs)
                
                # polygons_gdf = gpd.GeoDataFrame({'color': self.polygon_colors}, geometry=polygon_geometries, crs=self.left_rio.crs)
                polygons_gdf.to_file(filename, driver = 'GPKG')
        else:
            messagebox.showwarning("Warning", "No polygons drawn. Draw polygons first.")        

    def load_polygons(self):
        self.upload_polygons_flag = True
        self.best_denom_id_list = []
        self.ratio_spectrum_list = []
        # Function to parse and clean a parameter values vector string
        def parse_parameter_values(param_str):
            param_str = param_str.replace('[', '').replace(']', '').replace('\n', '')
            parameter_values = param_str.split()
            parameter_values = np.array([float(val) for val in parameter_values])
            return parameter_values
        
        # Function to parse and clean a spectrum string
        def parse_spectrum(spectrum_str):
            spectrum_str = spectrum_str.replace('[', '').replace(']', '').replace('\n', '')
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
                mean_spectrum = gpf['Spectrum Mean'].apply(parse_spectrum).tolist()
                parameter_values = gpf['Parameters Mean'].apply(parse_parameter_values).tolist()
                # if 'Ratio Spectrum' in gpf.columns:
                #     ratio_spectra = gpf['Ratio Spectrum'].apply(parse_spectrum).tolist()
            else:
                # Handle unsupported file types
                print("Unsupported file type. Please select a .shp or .gpkg file.")
                return
            self.restored_min_id_lines = []
            for i, geom in enumerate(gpf.geometry):
                print(f'reading in polygon {i}')
                # (lon, lat) <--> (x,y)
                polygon_coords = geom.exterior.coords.xy
                polygon_pixel_coords = []
                # translate from geographic to pixel coordinates
                for lon, lat in zip(polygon_coords[0], polygon_coords[1]):
                    y, x = self.left_rio.index(lon, lat)
                    polygon_pixel_coords.append((x, y))
                self.all_polygons.append(polygon_pixel_coords)
                self.polygon_colors.append(gpf.Color[i])
                self.restored_min_id_lines.append([gpf.get('Mineral ID 1')[i], gpf.get('Mineral ID 2')[i], gpf.get('Mineral ID 3')[i], gpf.get('Mineral ID 4')[i]])
                # if the gpf has 'Best Denom ID' column, populate the denominator index
                if 'Best Denom ID' in gpf.columns:
                    self.best_denom_id_list.append(gpf.get('Best Denom ID')[i])
                    print(f'Best Denom ID: {gpf.get("Best Denom ID")[i]}')
                else:
                    self.best_denom_id_list.append(np.nan)
                if 'Ratio Spectrum' in gpf.columns:
                    self.ratio_spectrum_list.append(gpf.get('Ratio Spectrum')[i])
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

        if hasattr(self, 'collected_points_tree'):
            if self.collected_points_window.winfo_exists():
                self.collected_points_tree.insert("","end", values = row)

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
        save_button = tk.Button(self.collected_points_window, text="Save to Excel", command=self.save_collected_points_to_excel)
        save_button.pack()

        # Add a button to plot the collected spectra
        plot_button = tk.Button(self.collected_points_window, text="Plot Spectra", command=self.plot_collected_spectra)
        plot_button.pack()

        self.collected_points_tree.bind("<Double-1>", self.edit_points_tree_item)
        self.collected_points_tree.bind("<Button-2>", self.show_collected_points_tree_context_menu)

        self.collected_points_tree_context_menu = tk.Menu(self.collected_points_tree, tearoff=0)
        self.collected_points_tree_context_menu.add_command(label="Delete point", command=self.delete_collected_point)
        
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
            self.entry.insert(0, self.collected_points_tree.item(selected_item, "values")[int(column[1:]) - 1])

            # Place the Entry widget over the cell
            self.entry.place(x=event.x, y=event.y, width=self.collected_points_tree.column(column, "width"))

            # Bind the Return key and the focus out event to the save_edit method
            self.entry.bind("<Return>", lambda _: self.save_collected_point_table_edit(selected_item, column))
            self.entry.bind("<FocusOut>", lambda _: self.save_collected_point_table_edit(selected_item, column))

    def save_collected_point_table_edit(self, item, column):
        # Get the new value
        new_value = self.entry.get()

        # Update the Treeview
        self.collected_points_tree.set(item, column, new_value)

        # Update the DataFrame
        self.collected_points.loc[self.collected_points["name"] == self.collected_points_tree.item(item, "values")[0], self.collected_points_tree.column(column, "id")] = new_value

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
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])

        # Save the DataFrame to an Excel file
        if file_path:
            self.collected_points.to_excel(file_path)

    def delete_collected_point(self):
        selected_item = self.collected_points_tree.selection()[0]  # get selected item
        item_values = self.collected_points_tree.item(selected_item, "values")
        point_name = item_values[0]  # assuming the first column is the point name
        self.collected_points = self.collected_points[self.collected_points.name != point_name]
        self.collected_points_tree.delete(selected_item)
        
    def collect_point(self):
        column_headers = self.collected_points.columns.tolist()
        point_index = len(self.collected_points)
        point_name = f'{self.img_base_name}_{point_index}'
        x, y = self.current_spec_loc
        spec = (self.left_wvl, self.spectrum)
        params = (self.right_wvl, self.parameter_values)
        # get geolocation from x, y
        lon, lat = self.left_rio.xy(y, x)
        if self.denom_spectrum is not None:
            denom = self.denom_spectrum
        else:
            denom = ""
        row = (point_name,"", "", "", "", "", x, y, lon, lat, spec, denom, params)
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
        self.collected_spectra_canvas = FigureCanvasTkAgg(fig, master=self.collected_spectral_plot_window)

        # Create an Axes object and plot the spectra
        ax = fig.add_subplot(111)
        for _, row in self.collected_points.iterrows():
            # retrieve the values from the "spectrum" column
            spec_label = row["name"]
            wvl, spectrum = row["spectrum"]
            color = row["color"]
            if color == "":
                color = 'k'
            ax.plot(wvl, spectrum, label=spec_label, color=color)

        # add a button to turn the legend on or off
        self.toggle_collected_legend_button = tk.Button(ui_frame, text="Legend (On)", command=self.toggle_collected_spectral_plot_legend)
        self.toggle_collected_legend_button.pack(side=tk.TOP)

        # Add a button to reset x-axis span
        self.reset_spectral_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_collected_x_axis_span)
        self.reset_spectral_x_axis_button.pack(side=tk.TOP)

        # Add built-in span options to a dropdown menu
        span_options = ["Full Span", "410 - 1000 nm", "410 - 2500 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
        self.collected_span_var = tk.StringVar()
        self.collected_span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(ui_frame, textvariable=self.collected_span_var, values=span_options, state="readonly")
        span_menu.pack(side=tk.TOP)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_collected_x_axis_span)

        # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        self.collected_library_spectra_var = tk.StringVar(ui_frame)
        self.collected_library_spectra_var.set("None") # default value
        self.collected_library_spectra_list = [name.split('/')[-1] for name in self.usgs_spectra_folders]
        # sort library_spectra_list alphabetically
        self.collected_library_spectra_list.sort()
        # label the library spectra drop down menu "library spectra"
        self.collected_library_spectra_label = tk.Label(ui_frame, text="Library Spectra")
        self.collected_library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.collected_library_spectra_dropdown = tk.OptionMenu(ui_frame, self.collected_library_spectra_var, *self.collected_library_spectra_list, command=self.plot_collected_library_spectra)
        self.collected_library_spectra_dropdown.pack(side=tk.TOP)

        # add a button to remove library spectra from the plot
        self.remove_collected_library_spectra_button = tk.Button(ui_frame, text="Remove Library Spectra", command=self.remove_collected_library_spectra)
        self.remove_collected_library_spectra_button.pack(side=tk.TOP)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(self.collected_spectra_canvas, self.collected_spectral_plot_window)
        toolbar.update()
        self.collected_spectra_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Add span selector for x-axis
        self.create_collected_x_axis_span_selector()

        # Pack the FigureCanvasTkAgg object into the window
        self.collected_spectra_canvas.draw()
        self.collected_spectra_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    def toggle_collected_spectral_plot_legend(self):
        ax = self.collected_spectra_canvas.figure.axes[0]
        if ax.get_legend() is None:
            ax.legend()
            self.toggle_collected_legend_button.config(text = "Legend (On)", relief="raised")
        else:
            ax.legend_.remove()
            self.toggle_collected_legend_button.config(text = "Legend (Off)", relief="sunken")
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
            ax, on_x_span_select, 'horizontal', useblit=True)

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
            "2700 - 3900 nm": (2700, 3900)
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
            path_to_spec_data = self.usgs_spectra_path + event + '/' + event + '.txt'
            path_to_wvl_data = self.usgs_spectra_path + event + '/' + '*Wavelengths*' + '.txt'
            library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
            library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            # plot the values in spectrum_df on self.spectral_ax 
            ax.plot(library_wvl, library_reflectance, label=self.collected_library_spectra_var.get())
            # scale each line so that the min and max values equal the min and max max value of self.spectrum
            for line in ax.lines[1:]:
                new_y_data = line.get_ydata() * (np.nanmax(self.spectrum) / np.nanmax(line.get_ydata()) - np.nanmin(self.spectrum))+ np.nanmin(self.spectrum)
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

            polygons_spectral_figure, self.polygons_spectral_ax = plt.subplots(figsize=(5,3))

            for poly_color, s in zip(self.polygon_colors, self.polygon_spectra):
                s = np.array(s)
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
            self.reset_polygons_x_axis_button.pack(side=tk.TOP)

            # Add a button to ratio spectra to a a chosen denominator spectrum
            self.ratio_spectra_button = tk.Button(ui_frame, text="Ratio Spectra", command=self.ratio_polygon_spectra)
            self.ratio_spectra_button.pack(side=tk.TOP)

            # Add built-in span options to a dropdown menu
            span_options = ["Full Span", "410 - 1000 nm", "410 - 2500 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
            self.polygons_span_var = tk.StringVar()
            self.polygons_span_var.set("Full Span")  # Set the default span option
            span_menu = ttk.Combobox(ui_frame, textvariable=self.polygons_span_var, values=span_options, state="readonly")
            span_menu.pack(side=tk.RIGHT)

            # add a button to turn the legend on or off
            self.toggle_polygons_legend_button = tk.Button(ui_frame, text="Legend (On)", command=self.toggle_polygons_spectral_plot_legend)
            self.toggle_polygons_legend_button.pack(side=tk.TOP)

            # add a button to ignore bad bands
            self.polygons_ignore_bad_bands_button = tk.Button(ui_frame, text="Ignore Bad Bands (Off)", command=self.toggle_polygons_ignore_bad_bands)
            self.polygons_ignore_bad_bands_button.pack(side=tk.TOP)

            # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
            self.polygons_library_spectra_var = tk.StringVar(ui_frame)
            self.polygons_library_spectra_var.set("None") # default value
            self.polygons_library_spectra_list = [name.split('/')[-1] for name in self.usgs_spectra_folders]
            # sort library_spectra_list alphabetically
            self.polygons_library_spectra_list.sort()
            # label the library spectra drop down menu "library spectra"
            self.polygons_library_spectra_label = tk.Label(ui_frame, text="Library Spectra")
            self.polygons_library_spectra_label.pack(side=tk.TOP)
            # create the drop down menu
            self.polygons_library_spectra_dropdown = tk.OptionMenu(ui_frame, self.polygons_library_spectra_var, *self.polygons_library_spectra_list, command=self.polygons_plot_library_spectra)
            self.polygons_library_spectra_dropdown.pack(side=tk.TOP)

            # add a button to remove library spectra from the plot
            self.polygons_remove_library_spectra_button = tk.Button(ui_frame, text="Remove Library Spectra", command=self.polygons_remove_library_spectra)
            self.polygons_remove_library_spectra_button.pack(side=tk.TOP)

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
                s = np.array(s)
                s = s.flatten()
                if self.ignore_polygons_bad_bands_flag:
                    s = np.where(np.array(self.bad_bands)==0, np.nan, s)
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
    
    def toggle_polygons_spectral_plot_legend(self):
        ax = self.polygons_spectral_canvas.figure.axes[0]
        if ax.get_legend() is None:
            ax.legend()
            self.toggle_polygons_legend_button.config(text = "Legend (On)", relief="raised")
        else:
            ax.legend_.remove()
            self.toggle_polygons_legend_button.config(text = "Legend (Off)", relief="sunken")
        self.polygons_spectral_canvas.draw()

    def toggle_polygons_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_polygons_bad_bands_flag:
                self.ignore_polygons_bad_bands_flag = False
                self.polygons_ignore_bad_bands_button.config(text = "Ignore Bad Bands (Off)", relief="sunken")
            else:
                self.ignore_polygons_bad_bands_flag = True
                self.polygons_ignore_bad_bands_button.config(text = "Ignore Bad Bands (On)", relief="raised")
            self.update_polygons_spectral_plot()

    def polygons_plot_library_spectra(self, event):
        if self.polygons_library_spectra_var.get() != "None":
            # plot the selected library spectrum
            ax = self.polygons_spectral_canvas.figure.axes[0]
            path_to_spec_data = self.usgs_spectra_path + self.polygons_library_spectra_var.get() + '/' + self.polygons_library_spectra_var.get() + '.txt'
            path_to_wvl_data = self.usgs_spectra_path + self.polygons_library_spectra_var.get() + '/' + '*Wavelengths*' + '.txt'
            library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
            library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            # plot the values in spectrum_df on self.spectral_ax 
            ax.plot(library_wvl, library_reflectance, label=self.polygons_library_spectra_var.get())
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
                new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * ((y_max - y_min) / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))) + y_min
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
            "410 - 2500 nm": (410, 2500),
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
        # Polygon ratios
        # ----------------------------------------------------------------
    
    # ----------------------------------------------------------------
    #  Polygon ratio spectra
    # ----------------------------------------------------------------
    def ratio_polygon_spectra(self):
        if self.denominator_index < 0 or self.denominator_index >= len(self.polygon_spectra):
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

            polygons_ratio_spectral_figure, self.polygons_ratio_spectral_ax = plt.subplots(figsize=(5,3))

            for i, (poly_color, s) in enumerate(zip(self.polygon_colors, self.polygon_ratio_spectra)):
                if i == self.denominator_index: continue
                s = s.flatten()
                l = self.polygons_ratio_spectral_ax.plot(self.left_wvl, s, color=poly_color)
                self.polygon_ratio_spectral_lines.append(l)
            # self.polygon_ratio_spectral_lines, = self.polygons_spectral_ax.plot(self.left_wvl, self.polygon_spectra)
            xmin, xmax = self.polygons_ratio_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            all_min_y = []
            all_max_y = []
            for s in self.polygon_ratio_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_ratio_spectral_canvas = FigureCanvasTkAgg(polygons_ratio_spectral_figure, master=self.polygons_ratio_spectral_window)
            self.polygons_ratio_spectral_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

            # add a button to turn the legend on or off
            self.toggle_polygons_ratio_legend_button = tk.Button(ui_frame, text="Legend (On)", command=self.toggle_polygons_ratio_spectral_plot_legend)
            self.toggle_polygons_ratio_legend_button.pack(side=tk.TOP)

            # Add a button to reset x-axis span
            self.reset_polygons_ratio_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_polygons_ratio_x_axis_span)
            self.reset_polygons_ratio_x_axis_button.pack(side=tk.TOP)

            # Add built-in span options to a dropdown menu
            span_options = ["Full Span", "410 - 1000 nm", "410 - 2500 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
            self.polygons_ratio_span_var = tk.StringVar()
            self.polygons_ratio_span_var.set("Full Span")  # Set the default span option
            span_menu = ttk.Combobox(ui_frame, textvariable=self.polygons_ratio_span_var, values=span_options, state="readonly")
            span_menu.pack(side=tk.TOP)

            # Bind an event to update the x-axis span when a span option is selected
            span_menu.bind("<<ComboboxSelected>>", self.update_polygons_ratio_x_axis_span)

            # add a button to ignore bad bands
            self.polygons_ratio_ignore_bad_bands_button = tk.Button(ui_frame, text="Ignore Bad Bands (Off)", command=self.toggle_polygons_ratio_ignore_bad_bands)
            self.polygons_ratio_ignore_bad_bands_button.pack(side=tk.TOP)

            # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
            self.polygons_ratio_library_spectra_var = tk.StringVar(ui_frame)
            self.polygons_ratio_library_spectra_var.set("None") # default value
            self.polygons_ratio_library_spectra_list = [name.split('/')[-1] for name in self.usgs_spectra_folders]
            # sort library_spectra_list alphabetically
            self.polygons_ratio_library_spectra_list.sort()
            # label the library spectra drop down menu "library spectra"
            self.polygons_ratio_library_spectra_label = tk.Label(ui_frame, text="Library Spectra")
            self.polygons_ratio_library_spectra_label.pack(side=tk.TOP)
            # create the drop down menu
            self.polygons_ratio_library_spectra_dropdown = tk.OptionMenu(ui_frame, self.polygons_ratio_library_spectra_var, *self.polygons_ratio_library_spectra_list, command=self.polygons_ratio_plot_library_spectra)
            self.polygons_ratio_library_spectra_dropdown.pack(side=tk.TOP)

            # add a button to remove library spectra from the plot
            self.polygons_ratio_remove_library_spectra_button = tk.Button(ui_frame, text="Remove Library Spectra", command=self.polygons_ratio_remove_library_spectra)
            self.polygons_ratio_remove_library_spectra_button.pack(side=tk.TOP)

            # Create a toolbar for the spectral plot
            toolbar = NavigationToolbar2Tk(self.polygons_ratio_spectral_canvas, self.polygons_ratio_spectral_window)
            toolbar.update()
            self.polygons_ratio_spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            
            # Add span selector for x-axis
            self.create_polygons_ratio_x_axis_span_selector(self.left_wvl)
      
    def toggle_polygons_ratio_spectral_plot_legend(self):
        ax = self.polygons_ratio_spectral_canvas.figure.axes[0]
        if ax.get_legend() is None:
            ax.legend()
            self.toggle_polygons_ratio_legend_button.config(text = "Legend (On)", relief="raised")
        else:
            ax.legend_.remove()
            self.toggle_polygons_ratio_legend_button.config(text = "Legend (Off)", relief="sunken")
        self.polygons_ratio_spectral_canvas.draw()

    def toggle_polygons_ratio_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_bad_bands_flag:
                self.ignore_bad_bands_flag = False
                self.polygons_ratio_ignore_bad_bands_button.config(text = "Ignore Bad Bands (Off)", relief="sunken")
            else:
                self.ignore_bad_bands_flag = True
                self.polygons_ratio_ignore_bad_bands_button.config(text = "Ignore Bad Bands (On)", relief="raised")
            self.update_polygons_ratio_spectral_plot()
    
    def polygons_ratio_plot_library_spectra(self, event):
        if self.polygons_ratio_library_spectra_var.get() != "None":
            # plot the selected library spectrum
            ax = self.polygons_ratio_spectral_canvas.figure.axes[0]
            path_to_spec_data = self.usgs_spectra_path + self.polygons_ratio_library_spectra_var.get() + '/' + self.polygons_ratio_library_spectra_var.get() + '.txt'
            path_to_wvl_data = self.usgs_spectra_path + self.polygons_ratio_library_spectra_var.get() + '/' + '*Wavelengths*' + '.txt'
            library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
            library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            # plot the values in spectrum_df on self.spectral_ax 
            ax.plot(library_wvl, library_reflectance, label=self.polygons_ratio_library_spectra_var.get())
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
                new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * ((y_max - y_min) / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))) + y_min
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
        if self.polygons_ratio_spectral_window is None or not self.polygons_ratio_spectral_window.winfo_exists():
            self.polygon_ratio_spectral_lines = []
            self.create_polygons_ratio_spectral_plot()
        elif self.polygon_ratio_spectra:
            # clear everything from the plot before plotting
            self.polygons_ratio_spectral_ax.clear()
            self.polygon_ratio_spectral_lines = []

            for i, (c, s) in enumerate(zip(self.polygon_colors, self.polygon_ratio_spectra)):
                if i == self.denominator_index: continue
                s = s.flatten()
                if self.ignore_bad_bands_flag:
                    s = np.where(np.array(self.bad_bands)==0, np.nan, s)
                l = self.polygons_ratio_spectral_ax.plot(self.left_wvl, s, color=c)
                self.polygon_ratio_spectral_lines.append(l)
            xmin, xmax = self.polygons_ratio_spectral_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))

            all_min_y = []
            all_max_y = []
            for s in self.polygon_ratio_spectra:
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
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
            if i == self.denominator_index: continue
            min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
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
            "2700 - 3900 nm": (2700, 3900)
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
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
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
                if i == self.denominator_index: continue
                min_y, max_y = np.nanmin(s[xmin_idx:xmax_idx]), np.nanmax(s[xmin_idx:xmax_idx])
                all_min_y.append(min_y)
                all_max_y.append(max_y)
            min_y, max_y = np.nanmin(all_min_y), np.nanmax(all_max_y)

            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.polygons_ratio_spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

            self.polygons_ratio_spectral_canvas.draw()

        self.x_polygons_ratio_span_selector = SpanSelector(
            self.polygons_ratio_spectral_ax, on_x_span_select, 'horizontal', useblit=True)
        
    # ----------------------------------------------------------------
    # points
    # ----------------------------------------------------------------
    def create_spectral_plot(self):
        sw = spectral_window
        sw.create_spectral_plot(app)

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
            self.toggle_legend_button.config(text = "Legend (On)", relief="raised")
        else:
            self.spectral_ax.legend_.remove()
            self.toggle_legend_button.config(text = "Legend (Off)", relief="sunken")
        self.spectral_canvas.draw()

    def plot_library_spectra(self, event):
        # plot the selected library spectrum
        # self.point_plot_library_indices
        path_to_spec_data = self.usgs_spectra_path + event + '/' + event + '.txt'
        path_to_wvl_data = self.usgs_spectra_path + event + '/' + '*Wavelengths*' + '.txt'
        print(path_to_spec_data)
        library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
        library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
        xmin, xmax = self.spectral_ax.get_xlim()
        ymin, ymax = self.spectral_ax.get_ylim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        # plot the values in spectrum_df on self.spectral_ax 
        self.spectral_ax.plot(library_wvl, library_reflectance, label=self.library_spectra_var.get())
        # scale each line so that the min and max values equal the min and max max value of self.spectrum
        for line in self.spectral_ax.lines[1:]:
            new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * ((np.nanmax(self.spectrum[xmin_idx:xmax_idx]) - np.nanmin(self.spectrum[xmin_idx:xmax_idx])) / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))) + np.nanmin(self.spectrum[xmin_idx:xmax_idx])
            line.set_ydata(new_y_data)
        # min_y, max_y = np.nanmin(library_reflectance), np.nanmax(library_reflectance)
        # buffer = (max_y - min_y) * 0.1 
        # new_y_lim = (np.nanmin((ymin, min_y - buffer)), np.nanmax((ymax, max_y + buffer)))
        self.spectral_ax.set_ylim((ymin, ymax))
        self.spectral_ax.set_xlim(xmin, xmax)
        self.spectral_ax.legend()
        self.spectral_canvas.draw()

    def plot_mica_library_spectra(self, event):
        # plot the selected library spectrum
        # self.point_plot_library_indices
        path_to_spec_data = self.mica_spectra_path + event + '.tab'
        print(path_to_spec_data)
        mica_df = pd.read_csv(path_to_spec_data, header=None)
        print(mica_df)
        library_wvl = mica_df[0].values
        library_wvl = [val*1000.0 for val in library_wvl]
        library_reflectance = mica_df[1].values
        library_reflectance = np.where(library_reflectance>6000, np.nan, library_reflectance)
        print(library_wvl, library_reflectance)
        xmin, xmax = self.spectral_ax.get_xlim()
        ymin, ymax = self.spectral_ax.get_ylim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        # plot the values in spectrum_df on self.spectral_ax 
        self.spectral_ax.plot(library_wvl, library_reflectance, label=self.mica_library_spectra_var.get())
        # scale each line so that the min and max values equal the min and max max value of self.spectrum
        for line in self.spectral_ax.lines[1:]:
            new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * ((np.nanmax(self.spectrum[xmin_idx:xmax_idx]) - np.nanmin(self.spectrum[xmin_idx:xmax_idx])) / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))) + np.nanmin(self.spectrum[xmin_idx:xmax_idx])
            line.set_ydata(new_y_data)
            print(new_y_data)
        # min_y, max_y = np.nanmin(library_reflectance), np.nanmax(library_reflectance)
        # buffer = (max_y - min_y) * 0.1 
        # new_y_lim = (np.nanmin((ymin, min_y - buffer)), np.nanmax((ymax, max_y + buffer)))
        self.spectral_ax.set_ylim((ymin, ymax))
        self.spectral_ax.set_xlim(xmin, xmax)
        self.spectral_ax.legend()
        self.spectral_canvas.draw()

    def remove_library_spectra(self):
        # Remove all lines except the first one (assuming it's the main line of the plot)
        for line in self.spectral_ax.lines[1:]:
            line.remove()

        # Update the plot
        self.update_spectral_plot()

    def toggle_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_bad_bands_flag:
                self.ignore_bad_bands_flag = False
                self.ignore_bad_bands_button.config(text = "Ignore Bad Bands (Off)", relief="sunken")
            else:
                self.ignore_bad_bands_flag = True
                self.ignore_bad_bands_button.config(text = "Ignore Bad Bands (On)", relief="raised")
            self.update_spectral_plot()
        else:
            messagebox.showwarning("Warning", "No bad bands set. Set bad bands first. (Processing > Select Bad Bands)")

    def update_spectral_plot(self):
        if self.spectral_window is None or not self.spectral_window.winfo_exists():
            print('creating spectral plot')
            self.create_spectral_plot()
        else:
            if self.ignore_bad_bands_flag:
                self.spectrum = np.where(np.array(self.bad_bands)==0, np.nan, self.spectrum)
            self.spectral_line.set_ydata(self.spectrum)
            self.spectral_line.set_label(self.spectrum_label)
            self.spectral_ax.legend()
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
            "410 - 2500 nm": (410, 2500),
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
    # ratio point plots
    # ----------------------------------------------------------------
    def create_ratio_spectral_plot(self):
        self.ratio_spectral_window = tk.Toplevel(self.root)
        self.ratio_spectral_window.title("Ratio Spectral Plot")
        
        # Create a frame to hold UI elements with a fixed size
        ratio_ui_frame = tk.Frame(self.ratio_spectral_window)
        ratio_ui_frame.pack(fill=tk.X)
        
        # ratio_spectral_figure, self.ratio_spectral_ax = plt.subplots(figsize=(5,3))
        self.ratio_spectral_figure, (self.ratio_spectral_ax1, self.ratio_spectral_ax2) = plt.subplots(2, 1, figsize=(5, 6))
        self.ratio_spectral_axes = (self.ratio_spectral_ax1, self.ratio_spectral_ax2)

        # numerator and denominator are plotted on the top
        self.numerator_spectral_line, = self.ratio_spectral_ax1.plot(self.left_wvl, self.spectrum, label="numerator")
        self.denominator_spectral_line, = self.ratio_spectral_ax1.plot(self.left_wvl, self.denom_spectrum, label="denominator")
        xmin, xmax = self.ratio_spectral_ax1.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin((np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]))), np.nanmax((np.nanmax(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx])))
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.ratio_spectral_ax1.set_ylim(min_y - buffer, max_y + buffer)

        # ratioed spectrum plots on the bottom
        self.ratio_spectral_line, = self.ratio_spectral_ax2.plot(self.left_wvl, self.ratio_spectrum, label='ratio')
        xmin, xmax = self.ratio_spectral_ax2.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.ratio_spectral_ax2.set_ylim(min_y - buffer, max_y + buffer)

        self.ratio_spectral_ax1.legend(loc='best')  # 'loc' can be adjusted to specify the legend position
        self.ratio_spectral_ax2.legend(loc='best')  # 'loc' can be adjusted to specify the legend position

        self.ratio_spectral_canvas = FigureCanvasTkAgg(self.ratio_spectral_figure, master=self.ratio_spectral_window)
        self.ratio_spectral_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # ---
        # add a button to turn the legend on or off
        self.toggle_ratio_legend_button = tk.Button(ratio_ui_frame, text="Legend (On)", command=self.toggle_ratio_spectral_plot_legend)
        self.toggle_ratio_legend_button.pack(side=tk.TOP)

        # add a button to ignore bad bands
        self.ignore_bad_bands_button = tk.Button(ratio_ui_frame, text="Ignore Bad Bands (Off)", command=self.toggle_ratio_ignore_bad_bands)
        self.ignore_bad_bands_button.pack(side=tk.TOP)

        # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        self.ratio_library_spectra_var = tk.StringVar(ratio_ui_frame)
        self.ratio_library_spectra_var.set("None") # default value
        self.ratio_library_spectra_list = [name.split('/')[-1] for name in self.usgs_spectra_folders]
        # sort library_spectra_list alphabetically
        self.ratio_library_spectra_list.sort()
        # label the library spectra drop down menu "library spectra"
        self.ratio_library_spectra_label = tk.Label(ratio_ui_frame, text="Library Spectra")
        self.ratio_library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.ratio_library_spectra_dropdown = tk.OptionMenu(ratio_ui_frame, self.ratio_library_spectra_var, *self.ratio_library_spectra_list, command=self.plot_ratio_library_spectra)
        self.ratio_library_spectra_dropdown.pack(side=tk.TOP)

        # add a button to remove library spectra from the plot
        self.remove_ratio_library_spectra_button = tk.Button(ratio_ui_frame, text="Remove Library Spectra", command=self.remove_ratio_library_spectra)
        self.remove_ratio_library_spectra_button.pack(side=tk.TOP)

        # add a button to collect the point and save it to a table
        self.collect_ratio_point_button = tk.Button(ratio_ui_frame, text="Collect Point", command=self.collect_point)
        self.collect_ratio_point_button.pack(side=tk.TOP)
        # ----

        # Add a button to reset x-axis span
        self.reset_ratio_x_axis_button = tk.Button(ratio_ui_frame, text="Reset X-Axis Span", command=self.reset_ratio_x_axis_span)
        self.reset_ratio_x_axis_button.pack(side=tk.RIGHT)

        # Add built-in span options to a dropdown menu
        span_options = ["Full Span", "410 - 1000 nm", "410 - 2500 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
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

            # numerator and denominator are plotted on the top
            self.numerator_spectral_line.set_ydata(self.spectrum)
            self.denominator_spectral_line.set_ydata(self.denom_spectrum)
            xmin, xmax = self.ratio_spectral_ax1.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin((np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]))), np.nanmax((np.nanmax(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx])))
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax1.set_ylim(min_y - buffer, max_y + buffer)

            # ratioed spectrum plots on the bottom
            self.ratio_spectral_line.set_ydata(self.ratio_spectrum)
            xmin, xmax = self.ratio_spectral_ax2.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax2.set_ylim(min_y - buffer, max_y + buffer)

            self.ratio_spectral_ax1.legend(loc='best')  # 'loc' can be adjusted to specify the legend position
            self.ratio_spectral_ax2.legend(loc='best')  # 'loc' can be adjusted to specify the legend position

            self.ratio_spectral_canvas.draw()

    def plot_ratio_library_spectra(self, event):
        # plot the selected library spectrum
        # self.point_plot_library_indices
        path_to_spec_data = self.usgs_spectra_path + event + '/' + event + '.txt'
        path_to_wvl_data = self.usgs_spectra_path + event + '/' + '*Wavelengths*' + '.txt'
        print(path_to_spec_data)
        library_wvl = getWavelengthFromUSGS(path_to_wvl_data)
        library_reflectance = getReflectanceFromUSGS(path_to_spec_data)
        for i, ax in enumerate(self.ratio_spectral_axes):
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            # plot the values in spectrum_df on ax 
            ax.plot(library_wvl, library_reflectance, label=self.ratio_library_spectra_var.get())
            # scale each line so that the min and max values equal the min and max max value of self.spectrum
            for line in ax.lines[1:]:
                if i == 0:
                    s = self.spectrum
                elif i == 1:
                    s = self.ratio_spectrum
                new_y_data = (line.get_ydata() - np.nanmin(line.get_ydata())) * ((np.nanmax(s[xmin_idx:xmax_idx]) - np.nanmin(s[xmin_idx:xmax_idx])) / (np.nanmax(line.get_ydata()) - np.nanmin(line.get_ydata()))) + np.nanmin(s[xmin_idx:xmax_idx])
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
                self.toggle_ratio_legend_button.config(text = "Legend (On)", relief="raised")
            else:
                ax.legend_.remove()
                self.toggle_ratio_legend_button.config(text = "Legend (Off)", relief="sunken")
            self.ratio_spectral_canvas.draw()

    def toggle_ratio_ignore_bad_bands(self):
        # first check if self.bad_bands is set
        if hasattr(self, "bad_bands"):
            if self.ignore_ratio_bad_bands_flag:
                self.ignore_ratio_bad_bands_flag = False
                self.ignore_ratio_bad_bands_flag.config(text = "Ignore Bad Bands (Off)", relief="sunken")
            else:
                self.ignore_ratio_bad_bands_flag = True
                self.ignore_ratio_bad_bands_flag.config(text = "Ignore Bad Bands (On)", relief="raised")
            self.update_spectral_plot()
        else:
            messagebox.showwarning("Warning", "No bad bands set. Set bad bands first. (Processing > Select Bad Bands)")

    def reset_ratio_x_axis_span(self):
        # Reset x-axis span to the default range
        self.ratio_spectral_ax1.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.ratio_spectral_ax1.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin((np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]))), np.nanmax((np.nanmax(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx])))
        buffer = (max_y - min_y) * 0.1
        self.ratio_spectral_ax1.set_ylim(min_y - buffer, max_y + buffer)

        self.ratio_spectral_ax2.set_xlim(self.left_wvl[0], self.left_wvl[-1])
        xmin, xmax = self.ratio_spectral_ax2.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])
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
            "2700 - 3900 nm": (2700, 3900)
        }
        if selected_span in span_ranges:
            xlim = span_ranges[selected_span]
            self.ratio_spectral_ax1.set_xlim(xlim[0], xlim[1])
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[0]))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xlim[1]))
            min_y, max_y = np.nanmin((np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]))), np.nanmax((np.nanmax(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx])))
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
            y_min, y_max = np.nanmin((np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmin(self.denom_spectrum[xmin_idx:xmax_idx]))), np.nanmax((np.nanmax(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.denom_spectrum[xmin_idx:xmax_idx])))
            buffer = (y_max - y_min) * 0.1
            self.ratio_spectral_ax1.set_ylim(y_min - buffer, y_max + buffer)

            y_min = np.nanmin(self.ratio_spectrum[xmin_idx:xmax_idx])
            y_max = np.nanmax(self.ratio_spectrum[xmin_idx:xmax_idx])

            buffer = (y_max - y_min) * 0.1  # Add a buffer to y-limits
            self.ratio_spectral_ax2.set_ylim(y_min - buffer, y_max + buffer)

            self.ratio_spectral_canvas.draw()

        self.ratio_top_x_span_selector = SpanSelector(
            self.ratio_spectral_ax1, on_x_span_select, 'horizontal', useblit=True)

        self.ratio_btm_x_span_selector = SpanSelector(
            self.ratio_spectral_ax2, on_x_span_select, 'horizontal', useblit=True)

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
        self.bbl_spectrum = self.left_data[int(s1/2,), int(s2/2,), :].flatten()

        # create the spectral plot
        bad_bands_figure, self.bad_bands_ax = plt.subplots(figsize=(5,3))
        self.bad_bands_line, = self.bad_bands_ax.plot(self.left_wvl, self.bbl_spectrum)
        self.bad_bands_ax.set_xlabel('Wavelength')
        self.bad_bands_ax.set_ylabel('Value')
        self.bad_bands_ax.set_title('Spectral Plot')
        xmin, xmax = self.bad_bands_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.bbl_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.bbl_spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.bad_bands_ax.set_ylim(min_y - buffer, max_y + buffer)

        # place the plot in the bad_bands_window
        self.bad_bands_canvas = FigureCanvasTkAgg(bad_bands_figure, master=self.bad_bands_window)
        self.bad_bands_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # initialize the bad_bands list as all 1's (i.e., good)
        if hasattr(self, 'bad_bands'):
            if self.bad_bands is not None:
                self.bad_bands_line.set_ydata(np.where(np.array(self.bad_bands)==0, np.nan, self.bbl_spectrum))
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
            self.bad_bands_line.set_ydata(np.where(np.array(self.bad_bands)==0, np.nan, self.bbl_spectrum))
            self.bad_bands_canvas.draw()
            # update the bad band values in the table
            for i, value in enumerate(self.bad_bands):
                self.bad_bands_table.set(i, column='#2', value=value)

        self.bad_bands_x_span_selector = SpanSelector(
            self.bad_bands_ax, on_x_span_select, 'horizontal', useblit=True)
        
        # add a way to show the bad bands in a table
        self.bad_bands_table = ttk.Treeview(self.bad_bands_window, columns=('wavelength', 'bad'))
        self.bad_bands_table.heading('#0', text='Band')
        self.bad_bands_table.heading('wavelength', text='Wavelength')
        self.bad_bands_table.heading('bad', text='Bad=0')
        self.bad_bands_table.pack(fill=tk.BOTH, expand=True)

        # populate the table with the bad bands
        for i, value in enumerate(self.bad_bands):
            self.bad_bands_table.insert(parent='', index='end', iid=i, text=str(i+1),
                                        values=(self.left_wvl[i], value))

        # add a way to change the bad bands
        def on_bad_band_change(event):
            selection = self.bad_bands_table.selection()
            if selection:  # Check if there is a selection
                item = selection[0]
                column = self.bad_bands_table.identify_column(event.x)
                if column == '#2':  # bad column
                    wvl, v = self.bad_bands_table.item(item)['values']
                    if v==1:
                        self.bad_bands_table.set(item, column='#2', value=0)
                    elif v==0:
                        self.bad_bands_table.set(item, column='#2', value=1)
            # update self.bad_bands to reflect changes in the table
            self.bad_bands = [int(self.bad_bands_table.item(i)['values'][1]) for i in range(len(self.bad_bands))]
            # update the spectral plot to reflect the change in bad bands list
            self.bad_bands_line.set_ydata(np.where(np.array(self.bad_bands)==0, np.nan, self.bbl_spectrum))
            self.bad_bands_canvas.draw()

        self.bad_bands_table.bind('<Double-1>', on_bad_band_change)

        # add an option to reset the bad bands list to all good bands
        self.plot_spectrum_from_main_view = tk.Button(self.bad_bands_window, text="Plot Spectrum from Main View", command=self.get_main_view_spec)
        self.plot_spectrum_from_main_view.pack(side=tk.BOTTOM)

        # add an option to reset the bad bands list to all good bands
        self.reset_bad_bands_button = tk.Button(self.bad_bands_window, text="Reset Bad Bands", command=self.reset_bad_bands)
        self.reset_bad_bands_button.pack(side=tk.BOTTOM)

        # add an option to linearly interpolate across the bad bands for the whole image cube
        self.interpolate_bad_bands_button = tk.Button(self.bad_bands_window, text="Interpolate Bad Bands", command=self.interpolate_bad_bands)
        self.interpolate_bad_bands_button.pack(side=tk.BOTTOM)

        # add an option to restore the original cube
        self.restore_original_cube = tk.Button(self.bad_bands_window, text="Restore Uninterpolated Cube", command=self.restore_original_cube)
        self.restore_original_cube.pack(side=tk.BOTTOM)

    def reset_bad_bands(self):
        self.bad_bands = [1] * len(self.left_wvl)
        # update the spectral plot to reflect the change in bad bands list
        self.bad_bands_line.set_ydata(np.where(np.array(self.bad_bands)==0, np.nan, self.bbl_spectrum))
        self.bad_bands_canvas.draw()
        # update the bad band values in the table
        for i, value in enumerate(self.bad_bands):
            self.bad_bands_table.set(i, column='#2', value=value)

    def get_main_view_spec(self):
        # get the spectrum from the main view and plot it in the bad_bands_window
        if not self.spectrum.tolist():
            messagebox.showwarning("Warning", "No spectrum selected, click anywhere on the main view to load a spectrum.")
        else:
            self.bbl_spectrum = self.spectrum
            self.bad_bands_line.set_ydata(self.spectrum)
            xmin, xmax = self.bad_bands_ax.get_xlim()
            xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
            xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
            min_y, max_y = np.nanmin(self.bbl_spectrum[xmin_idx:xmax_idx]), np.nanmax(self.bbl_spectrum[xmin_idx:xmax_idx])
            buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
            self.bad_bands_ax.set_ylim(min_y - buffer, max_y + buffer)
            self.bad_bands_canvas.draw()

    def interpolate_bad_bands(self):
        self.left_data_og = self.left_data #store the original cube in case we need it.
        interp_cube = np.array(self.ld)
        # set bad bands to nan
        for i, band in enumerate(self.bad_bands):
            if band==0:
                interp_cube[:,:,i] = np.nan*interp_cube[:,:,i]
            elif band==1:
                pass
        
        print(np.shape(interp_cube))

        ni = interpNaNs(interp_cube, self.left_wvl)
        ni.linearInterp()
        self.left_data = ni.data_cube

    def restore_original_cube(self):
        if hasattr(self, 'left_data_og'):
            self.left_data = self.left_data_og
        else:
            messagebox.showwarning("Warning", "No original cube found. Run interpolation first.")    

    def calculate_spectral_parameters(self):
        self.spectral_parameters_window = tk.Toplevel(self.root)
        self.spectral_parameters_window.title("Spectral Parameterization")

        # Checkbox variables
        self.crop_var = tk.IntVar()
        self.bbl_var = tk.IntVar()
        self.interpNans_var = tk.IntVar()
        self.denoise_var = tk.IntVar()

        # Checkboxes
        tk.Checkbutton(self.spectral_parameters_window, text="Crop", variable=self.crop_var).pack()
        # Tooltip(self.spectral_parameters_window, "Crop: Specify crop region, like [row0, row1, column0, column1], with row and column values in pixel coordinates (not lat/lon)")

        tk.Checkbutton(self.spectral_parameters_window, text="BBL", variable=self.bbl_var, command=self.activate_bbl_function).pack()
        # Tooltip(self.spectral_parameters_window, "BBL: Bad Bands List, 1=good, 0=bad. Use Processing>Select Bad Bands to select bad bands.")

        tk.Checkbutton(self.spectral_parameters_window, text="Interpolate NaNs", variable=self.interpNans_var).pack()
        # Tooltip(self.spectral_parameters_window, "Interpolate NaNs: Option to interpolate NaNs")

        tk.Checkbutton(self.spectral_parameters_window, text="Denoise", variable=self.denoise_var).pack()
        # Tooltip(self.spectral_parameters_window, "Denois: Option to denoise, this may take a long time. Recommend using with Interpolate NaNs")

        # run button
        tk.Button(self.spectral_parameters_window, text="Run", command=self.run_spectral_parameters).pack()
        # run button
        tk.Button(self.spectral_parameters_window, text="Initiate", command=self.initiate_spectral_parameters).pack()

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
        self.pc = cubeParamCalculator(bbl=bbl, interpNans=interpNans, crop=crop, denoise=denoise)
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
        remove_button = Button(edit_window, text="Remove Parameter", command=remove_parameter)
        remove_button.pack()

        # Button to close the edit window
        close_button = Button(edit_window, text="Close", command=edit_window.destroy)
        close_button.pack()

    def activate_bbl_function(self):
        if not hasattr(self, 'bad_bands'):
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
            with open(file_path, 'wb') as file:
                pickle.dump(state_dict, file)

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
            self.polygon_params = restored_instance['polygon_params']
            self.spectrum = restored_instance['spectrum']
            self.polygons_table_data = restored_instance['polyong_table_data']
            self.ratio_spectrum = restored_instance['ratio_spectrum']
            self.bad_bands = restored_instance['bad_bands']
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

    def create_spectral_plot(self, pc_index = None, d = None):
        self.num_idx = pc_index
        if pc_index is not None:
            if d is not None:
                self.spectrum = self.polygon_spectra[pc_index] / self.polygon_spectra[self.denominator_index]
            else:
                self.spectrum = self.polygon_spectra[pc_index]
            self.spectrum_label = self.polygons_table_data[pc_index][0]
        self.spectral_window = tk.Toplevel(self.root)
        self.spectral_window.title("Spectral Plot")
        
        # Create a frame to hold UI elements with a fixed size
        ui_frame = tk.Frame(self.spectral_window)
        ui_frame.pack(side=tk.RIGHT, fill=tk.BOTH)

        spectral_figure, self.spectral_ax = plt.subplots(figsize=(5,3))
        self.spectral_line, = self.spectral_ax.plot(self.left_wvl, self.spectrum, label=self.spectrum_label)
        self.spectral_ax.legend(loc='best')  # 'loc' can be adjusted to specify the legend position
        self.spectral_ax.set_xlabel('Wavelength')
        self.spectral_ax.set_ylabel('Value')
        self.spectral_ax.set_title('Spectral Plot')
        
        xmin, xmax = self.spectral_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(self.spectrum[xmin_idx:xmax_idx]), np.nanmax(self.spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.spectral_ax.set_ylim(min_y - buffer, max_y + buffer)

        self.spectral_canvas = FigureCanvasTkAgg(spectral_figure, master=self.spectral_window)
        self.spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # add a button to save the spectrum
        self.store_ratio_spectrum = tk.Button(ui_frame, text="Add Spectrum to Table", command=self.add_spectrum_to_table)
        self.store_ratio_spectrum.pack(side=tk.TOP)

        # add a button to turn the legend on or off
        self.toggle_legend_button = tk.Button(ui_frame, text="Legend (On)", command=self.toggle_spectral_plot_legend)
        self.toggle_legend_button.pack(side=tk.TOP)

        # Add a button to reset x-axis span
        self.reset_spectral_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_x_axis_span)
        self.reset_spectral_x_axis_button.pack(side=tk.TOP)

        # Add built-in span options to a dropdown menu
        span_options = ["Full Span", "410 - 1000 nm", "410 - 2500 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
        self.span_var = tk.StringVar()
        self.span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(ui_frame, textvariable=self.span_var, values=span_options, state="readonly")
        span_menu.pack(side=tk.TOP)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_x_axis_span)

        # add a button to ignore bad bands
        self.ignore_bad_bands_button = tk.Button(ui_frame, text="Ignore Bad Bands (Off)", command=self.toggle_ignore_bad_bands)
        self.ignore_bad_bands_button.pack(side=tk.TOP)

        # USGS LIBRAY SPECTRA
        # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        self.library_spectra_var = tk.StringVar(ui_frame)
        self.library_spectra_var.set("None") # default value
        self.library_spectra_list = [name.split('/')[-1] for name in self.usgs_spectra_folders]
        # sort library_spectra_list alphabetically
        self.library_spectra_list.sort()
        # label the librar spectra drop down menu "library spectra"
        self.library_spectra_label = tk.Label(ui_frame, text="USGS Library Spectra")
        self.library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.library_spectra_dropdown = tk.OptionMenu(ui_frame, self.library_spectra_var, *self.library_spectra_list, command=self.plot_library_spectra)
        self.library_spectra_dropdown.pack(side=tk.TOP)

        # MICA LIBRARY SPECTRA
        # add a drop down menu to plot the library spectra, contents of the drop down menu are the library spectra located in the self.mica_library_spectra folder
        self.mica_library_spectra_var = tk.StringVar(ui_frame)
        self.mica_library_spectra_var.set("None") #default value
        self.mica_library_spectra_list = [name.split(".")[0].split('/')[-1] for name in self.mica_spectra_names]
        # sort mica_library_spectra_list alphabetically
        self.mica_library_spectra_list.sort()
        self.mica_library_spectra_label = tk.Label(ui_frame, text="MICA Library Spectra")
        self.mica_library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.mica_library_spectra_dropdown = tk.OptionMenu(ui_frame, self.mica_library_spectra_var, *self.mica_library_spectra_list, command=self.plot_mica_library_spectra)
        self.mica_library_spectra_dropdown.pack(side=tk.TOP)

        # add a button to remove library spectra from the plot
        self.remove_library_spectra_button = tk.Button(ui_frame, text="Remove Library Spectra", command=self.remove_library_spectra)
        self.remove_library_spectra_button.pack(side=tk.TOP)

        # add a button to collect the point and save it to a table
        self.collect_point_button = tk.Button(ui_frame, text="Collect Point", command=self.collect_point)
        self.collect_point_button.pack(side=tk.TOP)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(self.spectral_canvas, self.spectral_window)
        toolbar.update()
        self.spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Add span selector for x-axis
        self.create_x_axis_span_selector(self.left_wvl)

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
        single_spectrum_figure, self.single_spectrum_ax = plt.subplots(figsize=(5,3))
        self.single_spectrum_line, = self.single_spectrum_ax.plot(self.left_wvl, spectrum, label=self.polygons_table_data[self.pc_index][0])
        self.single_spectrum_ax.legend(loc='best')

        self.single_spectrum_ax.set_xlabel('Wavelength')
        self.single_spectrum_ax.set_ylabel('Value')


        xmin, xmax = self.single_spectrum_ax.get_xlim()
        xmin_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmin))
        xmax_idx = np.argmin(np.abs(np.array(self.left_wvl) - xmax))
        min_y, max_y = np.nanmin(spectrum[xmin_idx:xmax_idx]), np.nanmax(spectrum[xmin_idx:xmax_idx])
        buffer = (max_y - min_y) * 0.1  # Add a buffer to y-limits
        self.single_spectrum_ax.set_ylim(min_y - buffer, max_y + buffer)

        self.single_spectral_canvas = FigureCanvasTkAgg(single_spectrum_figure, master=self.single_spectrum_window)
        self.single_spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # add a button to turn the legend on or off
        self.toggle_legend_button = tk.Button(ui_frame, text="Legend (On)", command=self.toggle_spectral_plot_legend)
        self.toggle_legend_button.pack(side=tk.TOP)

        # Add a button to reset x-axis span
        self.reset_spectral_x_axis_button = tk.Button(ui_frame, text="Reset X-Axis Span", command=self.reset_x_axis_span)
        self.reset_spectral_x_axis_button.pack(side=tk.TOP)

        # Add built-in span options to a dropdown menu
        span_options = ["Full Span", "410 - 1000 nm", "410 - 2500 nm", "1000 - 2600 nm", "1200 - 2000 nm", "1800 - 2500 nm", "2000 - 2500 nm", "2700 - 3900 nm"]
        self.span_var = tk.StringVar()
        self.span_var.set("Full Span")  # Set the default span option
        span_menu = ttk.Combobox(ui_frame, textvariable=self.span_var, values=span_options, state="readonly")
        span_menu.pack(side=tk.TOP)

        # Bind an event to update the x-axis span when a span option is selected
        span_menu.bind("<<ComboboxSelected>>", self.update_x_axis_span)

        # add a button to ignore bad bands
        self.ignore_bad_bands_button = tk.Button(ui_frame, text="Ignore Bad Bands (Off)", command=self.toggle_ignore_bad_bands)
        self.ignore_bad_bands_button.pack(side=tk.TOP)

        # USGS LIBRARY SPECTRA
        # add a drop down menu to plot library spectra, contents of the drop down menu are the library spectra located in the library_spectra folder
        self.library_spectra_var = tk.StringVar(ui_frame)
        self.library_spectra_var.set("None") # default value
        self.library_spectra_list = [name.split('/')[-1] for name in self.usgs_spectra_folders]
        # sort library_spectra_list alphabetically
        self.library_spectra_list.sort()
        # label the librar spectra drop down menu "library spectra"
        self.library_spectra_label = tk.Label(ui_frame, text="Library Spectra")
        self.library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.library_spectra_dropdown = tk.OptionMenu(ui_frame, self.library_spectra_var, *self.library_spectra_list, command=self.plot_library_spectra)
        self.library_spectra_dropdown.pack(side=tk.TOP)

        # MICA LIBRARY SPECTRA
        # add a drop down menu to plot the library spectra, contents of the drop down menu are the library spectra located in the self.mica_library_spectra folder
        self.mica_library_spectra_var = tk.StringVar(ui_frame)
        self.mica_library_spectra_var.set("None") #default value
        self.mica_library_spectra_list = [name.split(".")[0].split('/')[-1] for name in self.mica_spectra_names]
        # sort mica_library_spectra_list alphabetically
        self.mica_library_spectra_list.sort()
        self.mica_library_spectra_label = tk.Label(ui_frame, text="MICA Library Spectra")
        self.mica_library_spectra_label.pack(side=tk.TOP)
        # create the drop down menu
        self.mica_library_spectra_dropdown = tk.OptionMenu(ui_frame, self.mica_library_spectra_var, *self.mica_library_spectra_list, command=self.plot_mica_library_spectra)
        self.mica_library_spectra_dropdown.pack(side=tk.TOP)

        # add a button to remove library spectra from the plot
        self.remove_library_spectra_button = tk.Button(ui_frame, text="Remove Library Spectra", command=self.remove_library_spectra)
        self.remove_library_spectra_button.pack(side=tk.TOP)

        # Create a toolbar for the spectral plot
        toolbar = NavigationToolbar2Tk(self.single_spectral_canvas, self.single_spectrum_window)
        toolbar.update()
        self.single_spectral_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Add span selector for x-axis
        self.create_x_axis_span_selector(self.left_wvl)

if __name__ == "__main__":
    root = tk.Tk()
    app = SpectralCubeAnalysisTool(root)
    root.protocol("WM_DELETE_WINDOW", app.on_closing)
    root.mainloop()
