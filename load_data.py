import spectral
import tkinter as tk
from tkinter import filedialog

class HyperspectralData:
    def __init__(self, cube_path):
        self.cube_path = cube_path
        self.load_cube()
    
    def load_cube(self):
        self.cube = spectral.io.envi.open(self.cube_path)
    
    def get_data(self):
        return self.cube.load()
    
    def get_wavelengths(self):
        return self.cube.bands.centers

    def get_band_count(self):
        return int(self.cube.metadata['bands'])

def select_cube_file():
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    cube_path = filedialog.askopenfilename(title="Select Hyperspectral Cube File", filetypes=[("HDR Files", "*.hdr")])
    return cube_path

def main():
    cube_path = select_cube_file()

    if not cube_path:
        print("No cube file selected. Exiting.")
        return
    
    hyperspectral_data = HyperspectralData(cube_path)
    
    data = hyperspectral_data.get_data()
    wavelengths = hyperspectral_data.get_wavelengths()
    band_count = hyperspectral_data.get_band_count()
    
    print(f"Loaded hyperspectral cube with {band_count} bands.")
    print(f"Wavelengths: {wavelengths}")

if __name__ == "__main__":
    main()
