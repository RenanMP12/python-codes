# Libraries 
import numpy as np
import matplotlib.pyplot as plt

# coordinate matrix
# coord = [nº of node | X | Y] 
coord = np.array([[1,    0, 33.3],
                  [2, 13.2, 62.3],
                  [3, 39.3, 84.5],
                  [4, 22.2, 30.1],
                  [5, 49.9, 57.6],
                  [6, 78.8, 78.2],
                  [7, 39.3, 10.0],
                  [8, 59.7, 34.3],
                  [9, 73.9, 36.2],
                 [10, 69.8,  5.1],
                 [11, 28.0, 50.0],
                 [12, 33.3, 55.0],
                 [13, 45.0, 49.0],
                 [14, 35.0, 30.0]])

fig, ax = plt.subplots(figsize = (12, 7))
plt.scatter(coord[:,1], coord[:,2])

# incidence matrix
# inci=[nº of element| node 1 | node 2 | node 3 | node 4 ] 
inci = np.array([[1, 1, 4, 5, 2],
                 [2, 2, 5, 6, 3],
                 [3, 4, 7, 8, 5],
                 [4, 5, 8, 9, 6],
                [5, 7, 10, 9, 8],
             [6, 11, 14, 13, 12]])

# pluviometer matrix
u = np.array([4.62, 3.81, 4.76, 5.45, 4.90, 10.35, 4.96, 4.26, 18.36, 15.69])

# plot element
fig, ax = plt.subplots(figsize = (12, 7))
for i in range(np.size(inci,0)):
    x = [coord[inci[i,1] - 1,1], coord[inci[i,2] - 1,1], coord[inci[i,3] - 1,1], coord[inci[i,4] - 1,1], coord[inci[i,1] - 1,1]]
    y = [coord[inci[i,1] - 1,2], coord[inci[i,2] - 1,2], coord[inci[i,3] - 1,2], coord[inci[i,4] - 1,2], coord[inci[i,1] - 1,2]]
    ax.fill(x,y)
    ax.title.set_text('Title')

# total area by analytical geometry
A = 0 

for i in range(np.size(inci, 1)):
    x1 = coord[inci[i,1] - 1,1]
    x2 = coord[inci[i,2] - 1,1]
    x3 = coord[inci[i,3] - 1,1]
    x4 = coord[inci[i,4] - 1,1]
    
    y1 = coord[inci[i,1] - 1,2]
    y2 = coord[inci[i,2] - 1,2]
    y3 = coord[inci[i,3] - 1,2]
    y4 = coord[inci[i,4] - 1,2]
    
    a = 0.5*(x1*y2-x2*y1-x1*y4+x2*y3-x3*y2+x4*y1+x3*y4-x4*y3)
    print('Element ' + str(i) + ': ' + str(a) + ' km²')
    A += a
    
print('Total area ' + str(A) + ' km²')

# total - rain
Q = 0
for i in range(np.size(inci, 1)):
    x1 = coord[inci[i,1] - 1,1]
    x2 = coord[inci[i,2] - 1,1]
    x3 = coord[inci[i,3] - 1,1]
    x4 = coord[inci[i,4] - 1,1]
    
    y1 = coord[inci[i,1] - 1,2]
    y2 = coord[inci[i,2] - 1,2]
    y3 = coord[inci[i,3] - 1,2]
    y4 = coord[inci[i,4] - 1,2]
    
    u1 = u[inci[i,1] - 1]
    u2 = u[inci[i,2] - 1]
    u3 = u[inci[i,3] - 1]
    u4 = u[inci[i,4] - 1]
    
    q = (4*(u1/4 - u2/4 - u3/4 + u4/4)*((x1*y3)/8 - (x3*y1)/8 - (x1*y4)/8 - (x2*y3)/8 + (x3*y2)/8 + (x4*y1)/8 + (x2*y4)/8 - (x4*y2)/8))/3 + (4*(u1/4 + u2/4 - u3/4 - u4/4)*((x1*y2)/8 - (x2*y1)/8 - (x1*y3)/8 + (x3*y1)/8 + (x2*y4)/8 - (x4*y2)/8 - (x3*y4)/8 + (x4*y3)/8))/3 + 4*(u1/4 + u2/4 + u3/4 + u4/4)*((x1*y2)/8 - (x2*y1)/8 - (x1*y4)/8 + (x2*y3)/8 - (x3*y2)/8 + (x4*y1)/8 + (x3*y4)/8 - (x4*y3)/8);
    print('Element ' + str(i) + ': ' + str(q) + ' mm³')
    
    Q += q
    
print('Total rain ' + str(Q) + ' mm³')
