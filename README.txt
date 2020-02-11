Input Uncertainty Assessment for 3D Geologic Modeling of Fault Zones

Author: Ashton Krajnovich (akrajnov@mines.edu)

Utilizes Monte Carlo style simulation to assess the uncertainty about a set of predefined geologic modeling inputs used for modeling fault zones. 

Associated with publication (in progress): "Input-based uncertainty assessment for poorly constrained 3D geologic modeling of fault zones"

The primary script is InputUncertaintyQuantification.py, which can be run to recreate the simulated data from the publication, including exporting and visualizing the data and useful simulation quality assessment plots. 

This script imports functions from a series of .py files (FUNCTION_*.py), and data from two different .csv files. 

Other essential data is initialized by the user where specified. 

Please go through the import * listings to ensure that all of the required packages are installed properly. Of note is the requirement of having an R distribution installed with the Rfast package installed. 

Code is commented throughout for ease of use, designed to be run in the Spyder IDE. 

This work is intended for research purposes only. Open source tools have been utilized throughout, and the author wishes to acknowledge the copyright and license requirements of the tools used. Also, the author wishes this tool to be used primarily for research purposes. The MIT license has been added to the repository to ensure this - if any authors of open source tools used in this work wish to discuss the use of this license, please notify the author. 

Work in progress - please notify the author (Ashton Krajnovich, akrajnov@mines.edu) regarding perceived issues or inefficiencies. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.