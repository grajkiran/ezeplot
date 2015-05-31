![Ezeplot](http://i.imgur.com/VAgBgLU.png)

Ezeplot - Dynamical systems visualisation (v1.0)
===============================================

Ezeplot is a tool for visualising autonomous dynamical systems. It can be used as an educational aid for studying non-linear dynamics.

Main features
------------
* Plotting and animating trajectories of (1d, 2d, 3d) autonomous dynamical systems
* 2d, 3d and polar modes
* Plotting the nullclines, time-series and direction field
* Detecting fixed points based on user-input
* Plotting Poincare sections

Obtaining Ezeplot
-----------------

### Minimal package
The minimal package is very small (< 200KB) and consists of just the Ezeplot code. Download this package if you already have recent versions of python (2.7 or 3.4), matplotlib (>=1.3.1), numpy(>=1.8.1) and scipy(>=0.13.3) installed on your computer.
Download the minimal set of files and run.

[ezeplot-1.0-minimal.zip](https://www.dropbox.com/s/48u1268ug92wyzv/ezeplot-1.0-minimal.zip?dl=1)

### Fully self-contained standalone package
The standalone package is much larger (~ 40MB) and includes Ezeplot bundled
along with python and the necessary libraries. Download, unzip and run - no
installation required.
|Platform | System requirements | Download |
|---------|---------------------|--------- |
|Windows  | Windows XP and above|[ezeplot-1.0-win64.zip](https://www.dropbox.com/s/7jqz1y8m6g7bv4j/ezeplot-1.0-win64.zip?dl=1)|
|Mac OS X | OS X 10.9 and above |[ezeplot-1.0-osx64.zip](https://www.dropbox.com/s/8c3jzdeb1h95fh5/ezeplot-1.0-osx64.zip?dl=1)|
|Linux    | Kernel 3.0 and above|[ezeplot-1.0-linux64.zip](https://www.dropbox.com/s/yt6ppi72rpnztfx/ezeplot-1.0-linux64.zip?dl=1)|

Running Ezeplot
---------------
### Minimal package
For the minimal package, you need to open the file __main__.py with your system's python executable.
On Windows, if python is associated with .py files, simply double clicking __main__.py is sufficient.
On Mac OS X and linux, change to the extracted directory, and type:

    python __main__.py

### Standalone package
For the standalone package, extract the contents of the zip file and run
ezeplot.bat, ezeplot.command or ezeplot.sh for windows, mac and linux
respectively.

Screenshots
-----------
![Lorentz attractor](http://i.imgur.com/0FaEiv3.png)
![Lorentz attractor - Time series](http://i.imgur.com/A0JeHsQ.png)
![Lorentz attractor - Poincare section](http://i.imgur.com/udQUbNs.png)

Credits
-------
* [Prof. Chirag Kalelkar](https://sites.google.com/site/kalelkar/), Department of Mechanical Engineering, Indian Institute of Technology, Kharagpur, India.

For encouragement, support and bug-hunting.

* Mr. Aditya Bandyopadhyay, Research Scholar, Department of Mechanical Engineering, Indian Institute of Technology, Kharagpur, India.

For beta-testing.

License
-------
Ezeplot is distributed under the GNU General Public License. See the file LICENSE for more details.

Links
-----
[Ezeplot](http://grajkiran.github.io/ezeplot)
[Python](https://www.python.org/)
[Matplotlib](http://matplotlib.org/)
[Scipy](http://scipy.org)
