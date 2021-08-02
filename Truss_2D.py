# Libraries
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

# Input
E = 80e9 # youngs modulus
A = 1e-4 # section area

# Coordinate matrix
coord = np.array([[1, 0, 0],[2, 1.2,0.3],[3, 2.4,0],[4,1.2,0.8]])
nnos = np.size(coord,0)  # number of nodes

# Incidence matrix
inci = np.array([[1, 1, 2],[2,2,3],[3,3,4],[4,4,1],[5,2,4]])
nel = np.size(inci,0) # number of elements

# Boundary conditions
#   bc=[node | degree of freedom | value]
#
#   Degree of freedom 1 --> x
#   Degree of freedom 2 --> y

bc = np.array([[1,1,0],[1,2,0],[3,2,0]])

# Mask stiffness matrix
mask = np.zeros((2*nnos,2*nnos))
for i in range(0, np.size(bc,0)):
    if bc[i,1] == 1:
        mask[2*bc[i,0]-2,2*bc[i,0]-2] = 1
    elif bc[i,1] == 2:
        mask[2*bc[i,0]-1,2*bc[i,0]-1] = 1
mask = ma.masked_equal(mask, 1)
mask = ma.mask_rowcols(mask)
mask = (mask==False)

# Mask load vector
maskv = np.zeros(2*nnos)
for i in range(0, np.size(bc,0)):
    if bc[i,1] == 1:
        maskv[2*bc[i,0]-2] = 1
    elif bc[i,1] == 2:
        maskv[2*bc[i,0]-1] = 1
maskv = ma.masked_equal(maskv, 1)
maskv = (maskv==False)

#  Load vector
#   F = [node | degree of freedom | value]
#
#   Degree of freedom 1 --> Fx
#   Degree of freedom 2 --> Fy

load = np.array([[4,1,-30e3],[4,2,20e3]])

f = np.zeros((2*nnos,1)) # global force vector pre-allocation

for i in range(0,np.size(load,0)):
    if load[i,1] == 1:
        f[int(2*load[i,0]-2)] = load[i,2]
    elif load[i,1] == 2:
        f[int(2*load[i,0]-1)] = load[i,2]


# Global matrix assembly
kg = np.zeros((2*nnos,2*nnos)) # global stiffness matrix pre-allocation
for i in range(0,nel):
    node1 = inci[i,1] # first node element
    node2 = inci[i,2] # second node element
    
    x1 = coord[int(node1) - 1, 1]
    x2 = coord[int(node2) - 1, 1]

    y1 = coord[int(node1) - 1, 2]
    y2 = coord[int(node2) - 1, 2]

    l = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    if (x2-x1) ==0:
        if y2 > y1:
            theta = 2*np.arctan(1)
        else:
            theta =-2*np.arctan(1)
    else:
        theta = np.arctan((y2-y1)/(x2-x1))
    
    c = np.cos(theta)
    s = np.sin(theta)
    
    # local stiffness matrix
    ke = E*A/l*np.array([[c**2, c*s, -c**2, -c*s], [c*s, s**2, -c*s, -s**2], [-c**2, -c*s, c**2, c*s], [-c*s, -s**2, c*s, s**2]])


    # localization vector
    loc = [2*node1-2,2*node1-1,2*node2-2,2*node2-1]
    
    # global stiffness matrix
    kg[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] = kg[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] +  ke
    
kg_aux = kg[mask.data]
kg_aux = np.reshape(kg_aux, (2*nnos-np.size(bc,0), 2*nnos-np.size(bc,0)))
f_aux  = f[maskv.data]
maskv.data

# Displacement
u = np.zeros((2*nnos, 1))
u[maskv.data] = np.linalg.solve(kg_aux, f_aux)
u_ = np.reshape(u,(nnos,2))
displacement_x = u_[:,0]
displacement_y = u_[:,1]

# Post-processing
factor = 5
new_coord = coord[:,1:3] + factor*u_
fig, ax = plt.subplots(figsize = (12, 7))

# original truss
for i in range(0, np.size(inci,0)):
    x1 = coord[int(inci[i,1]-1), 1]
    x2 = coord[int(inci[i,2]-1), 1]
    
    y1 = coord[int(inci[i,1]-1), 2]
    y2 = coord[int(inci[i,2]-1), 2]
    
    plt.plot([x1,x2], [y1,y2], 'bo-')
    
    s  = "{}".format(i+1)
    plt.text((x1+x2)/2, (y1+y2)/2 + 0.05, s, fontsize = 15)

# deformed truss
for i in range(0, np.size(inci,0)):
    x1 = new_coord[int(inci[i,1]-1), 0]
    x2 = new_coord[int(inci[i,2]-1), 0]
    
    y1 = new_coord[int(inci[i,1]-1), 1]
    y2 = new_coord[int(inci[i,2]-1), 1]
    
    plt.plot([x1,x2], [y1,y2], 'r--')    
    
# load applied
for i in range(0, np.size(load,0)):
    node = load[i,0]
    X = coord[int(node - 1), 1]
    Y = coord[int(node - 1), 2]
    U = 0
    V = 0
    if load[i,1] == 1:
        if load[i,2] > 0:
            U = 1
        else:
            U = -1
    elif load[i,1] == 2:
        if load[i,2] > 0:
            V = 1
        else:
            V = -1
    ax.quiver(X, Y, U, V)
    
fig.savefig('fem1a.png',dpi=300)    # save figure (dpi = number of pixels)

plt.title('Treliça', fontsize = 20)
plt.xlabel('Deslocamento horizontal', fontsize = 15)
plt.ylabel('Deslocamento vertical', fontsize = 15)
    
plt.show()

# Stress
stress = np.zeros((nel, 1))

for i in range(0,nel):
    node1 = inci[i,1] # first node element
    node2 = inci[i,2] # second node element

    x1 = coord[int(node1) - 1, 1]
    x2 = coord[int(node2) - 1, 1]

    y1 = coord[int(node1) - 1, 2]
    y2 = coord[int(node2) - 1, 2]

    l = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    if (x2-x1) ==0:
        if y2 > y1:
            theta = 2*np.arctan(1)
        else:
            theta =-2*np.arctan(1)
    else:
        theta = np.arctan((y2-y1)/(x2-x1))


    c = np.cos(theta)
    s = np.sin(theta)

    # local stiffness matrix
    ke = E*A/l*np.array([[c**2, c*s, -c**2, -c*s], [c*s, s**2, -c*s, -s**2], [-c**2, -c*s, c**2, c*s], [-c*s, -s**2, c*s, s**2]])

    # localization vector
    loc = [2*node1-2, 2*node1-1, 2*node2-2, 2*node2-1]

    d = u[[int(loc[0]),int(loc[1]),int(loc[2]), int(loc[3])]]

    elforce = np.matmul(ke,d)

    stress[i] = np.sqrt(elforce[1]**2 + elforce[2]**2)/A
    
    force_element = np.sqrt(elforce[1]**2 + elforce[2]**2)
    
    print("Normal force in element {0:5d} = {1:8.2f} N.".format(i + 1, force_element[0]))

    if ((x2-x1)*elforce[3])<0:
            stress[i]=-stress[i]


# Reaction forces
F = np.matmul(kg,u)
maskv = (maskv==False)
print(F[maskv.mask])
