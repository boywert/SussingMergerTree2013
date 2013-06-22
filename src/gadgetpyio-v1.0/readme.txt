This a simple module to read in the positions and velocities from
gadget type 1 snapshots.  It reads everything in code units, so the
user is responsible for knowing what they are doing.  

Created by:
Matthew Becker, University of Chicago, 2009
Matthew Turk, somewhere, 2009
Peter Teuben, somewhere else, 2009

This uses some of the example code from the public Gadget package.

To install, do this:

python setup.py build
python setup.py install

Some notes on installation:
1) If you are a mac user, you need to use the sudo command for the
last step above (i.e. sudo python setup.py install).
2) You need to have numpy.  I am not sure about compatibility, but I
think any version will work.

Send comments or bugs to:

Matthew Becker, becker.mr@gmail.com
Matthew Turk, matthewturk@gmail.com
Peter Teuben, teuben@astro.umd.edu
