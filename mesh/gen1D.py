import numpy as np
import sys

x0 = float(sys.argv[1])
xf = float(sys.argv[2])
nelem = int(sys.argv[3])
'''
nummat = int(sys.argv[4])
if (nummat < 1):
    sys.exit('Invalid material number')
if (nummat >= 1):
    j = 5
    while (j < (nummat*3)):
'''

x = np.linspace(x0, xf, nelem)
y = np.zeros(nelem+1)
flagnode = np.zeros(nelem)
flagnode[0] = 1
flagnode[-1] = 2
mat = np.zeros(nelem) + 1

# Node file
node = open('case1.n','wt')
print >>node, nelem
for i in range(nelem):
    print >>node, str(i)+':', x[i], y[i], int(flagnode[i])
node.close()

# Element file
elem = open('case1.e','wt')
print >>elem, nelem-1
for i in range(nelem-1):
    if (i<nelem-2):
        nn = i + 1
    else:
        nn = -1
    print >>elem, str(i)+':', i, i+1, -1, i-1, nn, -1, -1, -1, -1, (x[i]+x[i+1])/2.0, 0.0, int(mat[i])
elem.close()
