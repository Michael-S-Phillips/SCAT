PDS_VERSION_ID         = PDS3 

RECORD_TYPE            = STREAM
OBJECT                 = TEXT
  PUBLICATION_DATE     = 2015-04-28
  NOTE                 = "Introduction to the CRISM Type Spectra Library."
END_OBJECT             = TEXT
END 
                                                

============================================================================
Introduction
============================================================================

This archive contains 31 ratioed type spectra from 28 different Mars 
Reconnaissance Orbiter (MRO) Compact Reconnaissance Imaging Spectrometer for 
Mars (CRISM) targeted hyperspectral observations.

The primary documentation for this archive is the data set description found
in the file DATASET.CAT in the CATALOG directory, and in the publication
Viviano-Beck, C.E., et al., Revised CRISM Spectral Parameters and Summary 
Products Based on the Currently Detected Mineral Diversity on Mars, J. 
Geophys. Res., Vol. 119, pp. 1403-1431, doi:10.1002/2014JE004627, 2014.

Additional references are listed in the file REF.CAT in the CATALOG directory.

============================================================================
File Formats
============================================================================

The data files on this volume are provided in ASCII table format 
(extension .TAB). Each data file is accompanied by a detached PDS label file 
(extension .LBL), an ASCII file containing information describing the data 
file and data record format.

Document and catalog files (extension .CAT) are provided in ASCII format.

Image files are provided in Portable Network Graphics (.PNG) format.

Text files and detached label files contain 80-byte fixed-length ASCII
records, with a carriage-return character (ASCII 13) in the 79th byte and a
line-feed character (ASCII 10) in the 80th byte. This format allows the
files to be read by MacOS, DOS, UNIX and VMS operating systems. Note that
while both label and catalog files are ASCII, some operating systems will
not recognize the .LBL and .CAT extensions as ASCII types. It may be
necessary to either remap these extensions as ASCII in the system or text
browser, or else to open the files from within the text browser.


============================================================================
Volume Contents
============================================================================

The volume contains the following directories:

BROWSE - The files in the Browse directory are PNG images of plots of the 
tables in the Data directory. Each image has a detached PDS label.

CATALOG - All files in the Catalog directory are ASCII text documentation
of the data, mission, instrument, personnel and references relevant to this 
archive. They are also used to load the PDS Catalog so that the information 
is available for online searching. The file VOLDESC.CAT in the root 
directory is a catalog file that describes the archive volume. This is also 
an ASCII text file.

DATA - All data products are found in the Data directory. The Type Spectra 
Library data are all ASCII table files with detached PDS labels.

DOCUMENT - The Document directory contains the documentation necessary for 
understanding the archive. Included are the keyword definitions and type 
spectra classes documents.

INDEX - The Index directory contains the file INDEX.TAB, a table with one 
row for each file in the Data directory. The table lists the path and file 
name, the product ID, the product creation time, and other information about 
the product. The file INDEX.LBL is a PDS label that describes the table. All 
files in this directory are ASCII text.


============================================================================
Contacts
============================================================================

For questions concerning these data products and documentation, contact:

          PDS Geosciences Node
          Washington University
          Dept. of Earth and Planetary Sciences
          1 Brookings Drive
          St. Louis, MO 63130
          314-935-5493

          WWW Site:  http://pds-geosciences.wustl.edu
          Electronic mail address:  geosci@wunder.wustl.edu


For questions regarding PDS Standards or other archive volumes available
from the PDS, please contact PDS Operator at the PDS Central Node (at JPL):


    Internet            pds_operator@jpl.nasa.gov
    Telephone           (818) 354-4321
    U.S. Mail           Planetary Data System, PDS Operator
                        Jet Propulsion Laboratory
                        Mail Stop 202-101
                        4800 Oak Grove Dr.
                        Pasadena, CA 91109-8099
