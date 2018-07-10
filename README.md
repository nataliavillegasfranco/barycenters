# Barycenters
# Python implementation of crime response planning by linear programming

Response times for police calls are one of the most important factors in regards to dealing with crime. These times could be reduced by placing the police in strategic locations. In this project, we apply spectral clustering to partition the locations of crime incidents in the Denver area into smaller parts. Then, we implement and compare three discrete barycenter models to determine an optimal police distribution in each part of the partition. A google maps html file displays the crime hotspots and suggested police placement for the whole region.
======================================

How to compile and run?
-----------------------

*Clone the repository with the files:

barycenters.py	
fixed_transport.mod	
fixed_transport.run	
gmap.pickle	
murder.csv 

*You need Python 3.6.2 and to install the packages:

-pandas
-os
-numpy
-gmplot
-pickle
-subprocess
-re
-datetime
-calendar
-time
-sys

*You need to install AMPL:
https://ampl.com/

*Run:

barycenters.py


