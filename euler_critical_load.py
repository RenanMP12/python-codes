# Libraries
import numpy as np
import numpy.ma as ma
from scipy import linalg
np.set_printoptions(precision=2)

nel    = 2
nnos   = nel + 1
L      = 1
alldof = np.linspace(1, 1,2*nnos) # all degrees of freedom
kf = np.zeros((2*nnos,2*nnos))    # global stiffness matrix pre-allocation
kg = np.zeros((2*nnos,2*nnos))    # global stiffness matrix pre-allocation
coord = np.zeros((nnos, 2))       # coordinate matrix pre-allocation
inci = np.zeros((nel, 5))         # incidence matrix pre-allocation

for i in range(0, nnos):
    coord[i,0] = i + 1           # node number
    coord[i,1] = i*L/(nnos-1)    # node position

for i in range(0, nel):   
    inci[i,0] = i + 1          # element number
    inci[i,1] = i + 1          # first node
    inci[i,2] = i + 2          # second node
    inci[i,3] = coord[i,1]     # first node coordinate
    inci[i,4] = coord[i+1,1]   # second node coordinate

# Material properties
E = 200e9

# Geometry properties
# first element
D1 = 20e-3
d1 = 18e-3
# second element
D2 = 20e-3
d2 = 16e-3

l = L/nel

geo = np.array([[D1, d1],[D2, d2]])

# Boundary conditions
#   bc=[node | degree of freedom | value]
#
#   Degree of freedom 1 --> y
#   Degree of freedom 2 --> oz

bc = np.array([[1,1,0],[1,2,0]])

mask = np.zeros((2*nnos,2*nnos))
for i in range(0, np.size(bc,0)):
    if bc[i,1] == 1:
        mask[2*(bc[i,0] - 1),2*(bc[i,0] - 1)] = 1
    elif bc[i,1] == 2:
        mask[2*(bc[i,0] - 1)+1,2*(bc[i,0] - 1) + 1] = 1
mask = ma.masked_equal(mask, 1)
mask = ma.mask_rowcols(mask)
mask = (mask==False)

for i in range(nel):
    node1 = inci[i,1] # first node element
    node2 = inci[i,2] # second node element
    D = geo[i,0]
    d = geo[i,1]
    inertia = np.pi/64*(D**4 - d**4)
    
    # local stiffness matrix
    kf_e = E*inertia/l**3*np.array([[  12, 6*l, -12, 6*l],
                                    [ 6*l,  4*l**2, -6*l,  2*l**2],
                                    [-12,-6*l,  12,-6*l],
                                    [ 6*l,  2*l**2, -6*l,  4*l**2]])
        
    # local geometric matrix
    kg_e = 1/30/l*np.array([[ 36, 3*l, -36, 3*l],
                            [ 3*l,4*l**2, -3*l, -l**2],
                            [-36,-3*l,  36,-3*l],
                            [ 3*l,-l**2, -3*l, 4*l**2]])
        
    # localization vector
    loc = [2*node1-2,2*node1-1,2*node2-2,2*node2-1]
    
    # global stiffness matrix 
    kf[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] = kf[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] +  kf_e
    kg[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] = kg[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] +  kg_e
    

kg_ = kg[mask.data]
kg_ = np.reshape(kg_, (2*nnos-np.size(bc,0), 2*nnos-np.size(bc,0)))

kf_ = kf[mask.data]
kf_ = np.reshape(kf_, (2*nnos-np.size(bc,0), 2*nnos-np.size(bc,0)))

w, v = linalg.eig(kf_, kg_)
Pcr = np.min(np.real(w))

print('Pcr = ' + str(format(Pcr, '.2f')) + ' N')
