# Geologic-Model-Input-Uncertainty-Characterization

Geologic-Model-Input-Uncertainty-Characterization is a Python script and set of Python and R functions used to generate Monte Carlo samples of four specified geologic modeling inputs used in 3D geologic modeling of fault zones. The script allows for user-defined uncertainty parameterizations for key aspects of each chosen geologic modeling input (Fault trace, Structural orientation, Vertical termination depth and Fault zone thickness). It is associated with the publication "Uncertainty assessment for 3D geologic modeling of fault zones based on geologic inputs and prior knowledge" in Solid Earth by Krajnovich et al., 2020 (https://se.copernicus.org/preprints/se-2020-21/).

![Image of Yaktocat](https://github.com/ajkran2/Geologic-Model-Input-Uncertainty-Characterization/blob/master/fig07.png)

## Requirements

The contents of the repository should be cloned to a local disk location and necessary file handles marked throughout the main script (InputUncertaintyQuantification.py) should be updated according to your file system. The file can then be run in a python environment (tested in Spyder), given that the necessary dependencies are installed properly. 

## Dependencies

The easiest method of installing dependencies will be by installing the Anaconda distribution for Python 3.7 (https://www.anaconda.com/products/individual). Once installed, use the Anaconda Navigator software to install a version of RStudio as well, which will install a base R installation. 

The following is the list of non-default Python libraries that need to be installed in the Python environment. All libraries should be installed using conda when available. Links are provided to published libraries. 
- pymc3 (https://anaconda.org/conda-forge/pymc3)
- rpy2 (https://anaconda.org/r/rpy2)
- mplstereonet (https://pypi.org/project/mplstereonet/)

The following R library must also be installed (using RStudio).
- Directional (https://cran.r-project.org/web/packages/Directional/index.html)

## Usage

The primary script is InputUncertaintyQuantification.py, which can be run to recreate the simulated data from the publication, including exporting and visualizing the data and useful simulation quality assessment plots.

This script imports functions from a series of .py files (FUNCTION_*.py), and data from two different .csv files.

Other essential data is initialized by the user where specified.

Code is commented throughout for ease of use, and was developed and tested in the Spyder IDE.

## Disclaimers
This work is intended for research purposes only. Open source tools have been utilized throughout, and the author wishes to acknowledge the copyright and license requirements of the tools used. Also, the author wishes this tool to be used primarily for research purposes. The MIT license has been added to the repository to ensure this - if any authors of open source tools used in this work wish to discuss the use of this license, please notify the author.

Work in progress - please notify the author (Ashton Krajnovich, akrajnov@mines.edu) regarding perceived issues or inefficiencies.

## License
[MIT](https://choosealicense.com/licenses/mit/)

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
