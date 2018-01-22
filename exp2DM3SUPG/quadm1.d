#-----------#
# Unit quad case #
#-----------#

#=========
| POINTS |
=========#
5 # number of points #

# Nodes which define the boundary #
0:   0   0  0.05    1
1:   1   0  0.05    1
2:   1   1  0.05    1
3:   0   1  0.05    1

# Material markers #
4: 0.5   0.5    0     1 

#===========
| SEGMENTS |
===========#
4 # Number of segments #

# Boundary segments #
0:   0   1   1
1:   1   2   1
2:   2   3   1
3:   3   0   1
