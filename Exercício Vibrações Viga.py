#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numpy.ma as ma
from scipy import linalg
np.set_printoptions(precision=3)


# In[2]:


nel    = 1
nnos   = nel + 1
L      = 1
alldof = np.linspace(1, 1,2*nnos) # all degrees of freedom
kf = np.zeros((2*nnos,2*nnos))    # global stiffness matrix pre-allocation
m  = np.zeros((2*nnos,2*nnos))    # global stiffness matrix pre-allocation
coord = np.zeros((nnos, 2))       # coordinate matrix pre-allocation
inci = np.zeros((nel, 5))         # incidence matrix pre-allocation


# In[3]:


for i in range(0, nnos):
    coord[i,0] = i + 1           # node number
    coord[i,1] = i*L/(nnos-1)    # node position


# In[4]:


for i in range(0, nel):   
    inci[i,0] = i + 1          # element number
    inci[i,1] = i + 1          # first node
    inci[i,2] = i + 2          # second node


# In[5]:


# Material properties
E  = 70e9
ro = 3e3 
b  = 0.03
h  = 0.05

A = b*h
inertia = b*h**3/12

l = L/nel

# Boundary conditions
#   bc=[node | degree of freedom | value]
#
#   Degree of freedom 1 --> y
#   Degree of freedom 2 --> oz

bc = np.array([[1,1,0],[1,2,0]])


# In[6]:


mask = np.zeros((2*nnos,2*nnos))
for i in range(0, np.size(bc,0)):
    if bc[i,1] == 1:
        mask[2*(bc[i,0] - 1),2*(bc[i,0] - 1)] = 1
    elif bc[i,1] == 2:
        mask[2*(bc[i,0] - 1)+1,2*(bc[i,0] - 1) + 1] = 1
mask = ma.masked_equal(mask, 1)
mask = ma.mask_rowcols(mask)
mask = (mask==False)


# # Massa consistente

# In[7]:


for i in range(nel):
    node1 = inci[i,1] # first node element
    node2 = inci[i,2] # second node element
    
    # local stiffness matrix
    kf_e = E*inertia/l**3*np.array([[  12,     6*l,  -12,     6*l],
                                    [ 6*l,  4*l**2, -6*l,  2*l**2],
                                    [ -12,    -6*l,   12,    -6*l],
                                    [ 6*l,  2*l**2, -6*l,  4*l**2]])
    
    # local mass matrix
    m_e = ro*A*l/420*np.array([[   156,    22*l,    54,   -13*l],
                               [  22*l,  4*l**2,  13*l, -3*l**2],
                               [    54,    13*l,   156,   -22*l],
                               [ -13*l, -3*l**2, -22*l, 4*l**2]])
    
    # localization vector
    loc = [2*node1-2,2*node1-1,2*node2-2,2*node2-1]
    
    # global stiffness matrix 
    kf[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] = kf[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] +  kf_e
    m[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] = m[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] +  m_e
    
m_ = m[mask.data]
m_ = np.reshape(m_, (2*nnos-np.size(bc,0), 2*nnos-np.size(bc,0)))

kf_ = kf[mask.data]
kf_ = np.reshape(kf_, (2*nnos-np.size(bc,0), 2*nnos-np.size(bc,0)))


# In[8]:


w, v = linalg.eig(kf_, m_)


# # 1 modo

# In[9]:


w = np.min(np.real(w))
omega = np.sqrt(w)
f = omega/2/np.pi
print('omega = ' + str(format(omega, '.4f')) + ' rad/s')
print('f = ' + str(format(f, '.4f')) + ' Hz')


# In[10]:


(1.875)**2*np.sqrt(E*inertia/ro/A)


# # Massa concentrada

# In[11]:


kf = np.zeros((2*nnos,2*nnos))    # global stiffness matrix pre-allocation
m  = np.zeros((2*nnos,2*nnos))    # global stiffness matrix pre-allocation


# In[12]:


for i in range(nel):
    node1 = inci[i,1] # first node element
    node2 = inci[i,2] # second node element
    
    # local stiffness matrix
    kf_e = E*inertia/l**3*np.array([[  12,     6*l,  -12,     6*l],
                                    [ 6*l,  4*l**2, -6*l,  2*l**2],
                                    [ -12,    -6*l,   12,    -6*l],
                                    [ 6*l,  2*l**2, -6*l,  4*l**2]])
    
    # local mass matrix
    m_e = ro*A*l/2*np.array([[ 1, 0, 0, 0],
                             [ 0, 0, 0, 0],
                             [ 0, 0, 1, 0],
                             [ 0, 0, 0, 0]])
    
    # localization vector
    loc = [2*node1-2,2*node1-1,2*node2-2,2*node2-1]
    
    # global stiffness matrix 
    kf[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] = kf[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] +  kf_e
    m[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] = m[[[int(loc[0])],[int(loc[1])],[int(loc[2])], [int(loc[3])]], [int(loc[0]),int(loc[1]),int(loc[2]),int(loc[3])]] +  m_e
    
m_ = m[mask.data]
m_ = np.reshape(m_, (2*nnos-np.size(bc,0), 2*nnos-np.size(bc,0)))

kf_ = kf[mask.data]
kf_ = np.reshape(kf_, (2*nnos-np.size(bc,0), 2*nnos-np.size(bc,0)))


# In[13]:


w, v = linalg.eig(kf_, m_)


# In[14]:


w = np.min(np.real(w))
omega = np.sqrt(w)
f = omega/2/np.pi
print('omega = ' + str(format(omega, '.4f')) + ' rad/s')
print('f = ' + str(format(f, '.4f')) + ' Hz')

