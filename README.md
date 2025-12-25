# PyMAVEN
This package implements simple Python-based data analysis tools for NASA's [Mars Atmosphere & Volatile Evolution Mission (MAVEN)](https://pds-atmospheres.nmsu.edu/data_and_services) mission. This includes viewing data for the following MAVEN instruments: the Magnetometer (MAG), Solar Wind Electron Analyzer (SWEA), Solar Wind Ion Analyzer (SWIA), and SupraThermal and Thermal Ion Composition (STATIC). The interface uses [pyspedas](https://pyspedas.readthedocs.io/en/latest/) and [pytplot](https://pytplot.readthedocs.io/en/latest/) to load MAVEN data, but relies upon [matplotlib](https://matplotlib.org/) as its plotting backend for ease of customisation.

## Usage
To clone the git repo, run in a terminal window:
```sh
git clone https://github.com/kohzewen/pymaven.git
```
Then navigate to the new directory (you can use `cd pymaven`). To create a conda environment with the same dependencies used in this package, make sure you have `conda` installed and run:
```sh
conda env create -f pymaven.yaml
```
and
```sh
conda activate pymaven
```
You should now be able to run the example code `example.ipynb` in this package, which brings you through a few example functions you can use in this data anaysis package.

## Package contents
- `load_data.py`: Loads data from MAG, SWEA, SWIA, and STATIC instruments for a given date of interest.
- `process_data.py`:  Miscellaneous functions for data processing used in the backend.
- `plot_orbit.py`: Plots the orbit of MAVEN on a given day.
- `plot_inst.py`: Plots instrumental (MAG, SWEA, SWIA, STATIC) data in a given time interval.

## How to cite this package
If you use the PyMAVEN package, you can cite [this paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025GL117836), or use the following citation:
> Koh, Z., Poh, G., Fowler, C.M., Hanley, K.G., Ma, X., Gruesbeck, J. R., Kuruppuaratchi, D. C. P., Sun, W., DiBraccio, G. A., Espley, J. R., (2025). Global Occurrence of Kelvin-Helmholtz Vortices at Mars. Geophysical Research Letters, 52, e2025GL117836. DOI: 10.1029/2025GL117836

## License
[MIT](LICENSE) © Ze-Wen Koh
