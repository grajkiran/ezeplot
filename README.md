![Ezeplot](http://i.imgur.com/VAgBgLU.png)

Ezeplot - Dynamical systems visualisation (v1.0)
===============================================

Ezeplot is a tool for visualising autonomous dynamical systems. It can be used as an educational aid for studying non-linear dynamics.

Features
------------
* Plotting and animating trajectories of (1d, 2d, 3d) autonomous dynamical systems
* 2d, 3d and polar modes
* Plotting the nullclines, time-series and direction field
* Detecting fixed points based on user-input
* Plotting Poincare sections

Obtaining and running Ezeplot
-----------------------------

### Minimal package
The minimal package is very small (< 200KB) and consists of just the Ezeplot
code as a set of [python](https://www.python.org) programs.
Download this package if you already have recent versions of the following
software already installed on your computer:
* python (2.7 or 3.4)
* matplotlib (>=1.3.1)
* numpy(>=1.8.1)
* scipy(>=0.13.3)

You need to open the file `__main__.py` with your system's python executable.
On Windows, if python is associated with `*.py` files, opening `__main__.py` should work.

On Mac OSX and linux, change to the extracted directory, and type (in the terminal window):

    python __main__.py

Download the minimal package here: [ezeplot-1.0-minimal.zip](https://www.dropbox.com/s/48u1268ug92wyzv/ezeplot-1.0-minimal.zip?dl=1)

### Fully self-contained standalone package
The standalone package is much larger (~ 40MB) and includes Ezeplot bundled
along with python and the necessary libraries. This package requires no
installation and can also be run from a portable USB flash drive. 

#### Windows (Windows XP or above)
Download [ezeplot-windows.zip](https://www.dropbox.com/s/3ztkqm9xzhtzmw4/ezeplot-windows.zip?dl=1),
extract and run ezeplot.bat

#### Mac OSX (OSX 10.9 or above)
Download [ezeplot-osx.zip](https://www.dropbox.com/s/kmfstauqalyr5xb/ezeplot-osx.zip?dl=1),
extract and run ezeplot.command

#### Linux (Kernel 3.0 or above)
Download [ezeplot-linux.zip](https://www.dropbox.com/s/rbdtg0dayz6iacm/ezeplot-linux.zip?dl=1),
extract and run ezeplot.sh

Screenshots
-----------
![Lorentz attractor](http://i.imgur.com/Il9b2sf.png)
![Lorentz attractor - Time series](http://i.imgur.com/zAUQFlY.png)
![Lorentz attractor - Poincare section](http://i.imgur.com/BMeTMp2.png)

Credits
-------
* [Prof. Chirag Kalelkar](https://sites.google.com/site/kalelkar/), Department of Mechanical Engineering, Indian Institute of Technology, Kharagpur, India.

For encouragement, support and bug-hunting.

* Mr. Aditya Bandyopadhyay, Research Scholar, Department of Mechanical Engineering, Indian Institute of Technology, Kharagpur, India.

For beta-testing.

License
-------
Ezeplot is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl.html). See the file LICENSE for more details.

Links
-----
[Ezeplot](http://grajkiran.github.io/ezeplot)

[Python](https://www.python.org/)

[Matplotlib](http://matplotlib.org/)

[Scipy](http://scipy.org)
