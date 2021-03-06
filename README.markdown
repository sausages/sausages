# Sausages #

The `sausages` code efficiently analyses points on a regular three-dimensional grid to detect, analyse and classify continuous regions of adjacent points with similar properties.
Its original purpose was to detect and classify defect lines in liquid crystals, see e.g. [this paper](http://dx.doi.org/10.1080/02678292.2017.1295478) ([arXiv preprint](http://arxiv.org/abs/1702.02851)).
However, it could easily be expanded to apply to other systems.

This code was co-authored by Sam Brown (sausages@sambrown.eu) and Anja Humpert (anja.humpert@gmail.com).
It is released under the MIT license, and includes certain libraries with permissive licenses ([Eigen](http://eigen.tuxfamily.org/) for linear-algebra, [miniz](https://code.google.com/p/miniz/) for compression, and [cJSON](https://github.com/DaveGamble/cJSON) for JSON parsing).

## The physics ##
The figures below show two spherical colloidal particles surrounded by a liquid crystal.
The crystal is only displayed where it is in a particular (heavily biaxial) state, shown by the red ‘sausages’.
Biaxial grid-points correspond to a low value of a particular field (Cl), represented in the input file.
The biaxial regions "entangle" the colloids, and the different classes of entanglement are associated with different free energies (measured separately).

![Twin Rings defect structure](images/twin-rings.png)
![Omega defect structure](images/omega.png)
![Figure-of-Eight defect structure](images/figure-eight.png)
![Second Loop defect structure](images/second-loop.png)

The code firstly detects the red regions.
If the sausages can physically exist only as loops (as for liquid crystal defect lines), it then checks that all regions are closed loops;
if there are small gaps in the sausage due to noise, the gap is closed.
Finally the systems are classified into one of the above expected configurations.

## Compilation and execution ##
Clone the code from the Github repository:
```
git clone https://github.com/sausages/sausages.git
```
Some of the files used for end-to-end testing are quite large, and may cause time-out issues when cloning.
If this happens, try performing an initial shallow clone:
```
git clone --depth 1 https://github.com/sausages/sausages.git
git clone https://github.com/sausages/sausages.git
```

The code can be compiled using the provided Makefile by running:
```
make
```
The code includes C++11 elements that may create problems with older compilers; the code has been tested with GCC version 4.9.


To run an example file with example parameters try one of the tests:
```
cd tests/test1
../../sausages short.xyzclcpcs test.params
```

to run all the tests (advised):
```
cd tests
./run_all.sh
```

To compile the Doxygen documentation, run:
```
make docs
```
The default is to compile both HTML and LaTeX documentation.
If LaTeX is installed then the documentation can be compiled to PDF format by
```
make manual.pdf
```


## The structure of the program ##
* The input file is read in and all relevant data points are stored.
* The flood-fill algorithm is used to identify contiguous regions (sausages).
* The sausages are checked to determine whether they form closed loops.
* If they do not form closed loops the endpoints of each sausage are joined, subject to them being in close proximity without a third point nearby.
* The length of the loops is estimated.
* The structure of the loops are classified into one of the four different structures shown above, if possible.

### The input file ###
The program reads in files in the format shown below. Instead of providing the xyz data for each data point, it only specifies the system’s lower boundaries, the number of points and the spacing between points. Each point’s values are then read in column-major order.

```
> head -n 11 westin.diot
version 0.2
# total numVoxels: 5600000
# Assuming cubic
numVoxels 100 120 120
voxelSize 0.5 0.5 0.5
lowBounds -25 -30 -30
colloidPos 1 -25 0 0 0 # First number is index
beginClCpCs zyxInc
0.760129 0.024124 0.215748
0.760329 0.030349 0.209322
0.759563 0.037460 0.202977
```

In our example we focus solely on the property of the first column (Cl).
Only data points with Cl < threshold are read in and stored, other points are discarded to save memory.
The input files can be quite large and so the program automatically unzips files with the extension `.zip`.
In the `read_diot()` function the direct neighbours of each point are identified and can be accessed with pointers to any neighbouring particles.

### Identify contiguous regions using flood-fill ###
In the next step we use the flood-fill algorithm to find neighbouring points.
The `sausage` class corresponds to distinct (separate, non-contiguous) groups of contiguous points, together with connectivity information in the form of an irregular 3D multiply linked list.
In our calculation we ignore sausages that contain less than 1% of the total data points.
This minimum, as well as many other values, can be specified by the user in `test.params`, otherwise the default is used.


### Checking that sausages form closed loops ###
In liquid crystals defect sausages can only exist as loops across periodic boundaries.
For each sausage the program checks whether it forms a closed loop.
If they are already closed the next step is skipped.
Otherwise the program calculates the endpoints of all unclosed sausages by calculating the longest shortest-path distance between all points within the sausage.
Endpoints that are in close proximity with no third endpoint nearby are joined, and their sausages are merged.
The test for closed loops is repeated and only if all sausages are closed does the program proceed.

### Distinguishing twisted defect structures ###
In order to distinguish the different defect structures shown in the four pictures above, a flood-fill algorithm is used to determine the connectivity of four regions.
For each colloid we find a plane perpendicular to the line connecting the two colloids, centred on the colloid.
Points within a certain distance of these planes are added to regions depending on their location, seen below.

![Four regions for flood-fill based structure classification](images/four-regions.png)

Once again the flood-fill algorithm - restricted to move only towards the centre - is used to see which of these regions are connected.
Depending on the connections between regions the different structures can be successfully distinguished.

### Estimating the length of the defect line ###
Two methods are available to estimate the length of the sausage.

The first uses the all-pairs shortest-path Floyd-Warshall algorithm to find the two points separated by the longest path.
Assuming a well-formed sausage, the path separating these two points is half the sausage length.

The second is a ballistic "sphere-tracking" approach.
To move along the sausage we use a sphere that slowly increases in size until it contains 10 data points, whose position is averaged.
A new point close the sphere along the sausage is chosen and the process is repeated until the initial starting point is reached, with care never to reverse direction.
The distances between neighbouring averaged points is summed to estimate the total length of the sausage.
The red points in the figure below show the averaged points that were used for the length calculation.
The blue data points are the original data.

![Sphere-tracking length estimation](images/sphere-tracking.png)

### Classification ###
Finally, based on both the number of defect structures, their twist and their relative length, the four classes of entanglement can be distinguished.
Additionally, the handedness of the figure-of-eight twist can be found by examining projection onto the plane-of-best-fit of the loop.
