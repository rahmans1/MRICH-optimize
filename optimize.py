from skopt import gp_minimize
from skopt.plots import plot_convergence
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import linalg as lin
import sys

n=110
a= np.sqrt(np.pi*(950**2-100**2)/n)   # calculate cross-sectional size
b= 200                                    # 200 mm= zlength of each module

def Transform(x,y,z,a,b,c):  
  xyz=np.zeros(3)
  t=-1*np.arctan(y/z)
  p=np.arctan(x/z)
  xyz[0]= a*np.cos(p)+np.sin(p)*(b*np.sin(t)+c*np.cos(t))
  xyz[1]= b*np.cos(t)-c*np.sin(t)
  xyz[2]= -a*np.sin(p)+np.cos(p)*(b*np.sin(t)+c*np.cos(t))
  return [xyz]
  
  

def RotateOrigin(x,y,z,v, coor, a, b):   # Enter coordinate of center, vertex number, type of coordinate and half lengths along x,y,z
  xyz=None
  
  if v==0:
    xyz = Transform(x,y,z,-a,-a,-b)
  elif v==1:
    xyz = Transform(x,y,z,-a, a, -b)
  elif v==2:
    xyz = Transform(x,y,z, a, a, -b)
  elif v==3:
    xyz = Transform(x,y,z, a, -a, -b)
  elif v==4:
    xyz = Transform(x,y,z,-a, -a, b)
  elif v==5:
    xyz = Transform(x,y,z,-a, a, b)
  elif v==6:
    xyz = Transform(x,y,z, a, a, b)
  elif v==7:
    xyz = Transform(x,y,z, a, -a, b)

  if coor=='x':
    return xyz[0][0]
  elif coor=='y':
    return xyz[0][1]
  elif coor=='z':
    return xyz[0][2]
  else:
    return xyz[0]


def Normal(R):
   A1 = np.cross(R[7]-R[2], R[3]-R[2])
   A1 /= np.linalg.norm(A1)
   A2 = np.cross(R[6]-R[1], R[2]-R[1])
   A2 /= np.linalg.norm(A2)
   A3 = np.cross(R[5]-R[6], R[7]-R[6])
   A3 /= np.linalg.norm(A3)

 #  print(R)
 #  print(A1)
 #  print(A2)
 #  print(A3)

 #  print(np.dot(A1,A3))
 #  print(np.dot(A1,A2))
 #  print(np.dot(A2,A3))
   
   return [A1,A2,A3]


def Intersects(r1, r2, a=a/2 , b=b/2):    # Enter the position of the center for the two modules and half lengths along x,y, and z
   d = r2-r1   # Vector connecting center of two modules
   R1=[(r1+RotateOrigin(r1[0],r1[1],r1[2],v, 'all', a, b)) for v in range(0,8)]
   R2=[(r2+RotateOrigin(r2[0],r2[1],r2[2],v, 'all', a, b)) for v in range(0,8)]

   fig = go.Figure(data=[
      go.Mesh3d(
        # 8 vertices of a cube
        x=[R[v][0] for v in range(0,8)],
        y=[R[v][1] for v in range(0,8)],
        z=[R[v][2] for v in range(0,8)],
        colorbar_title='z',
        # i, j and k give the vertices of triangles
        i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
        j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
        name='y',
        showscale=True
      ) for R in [R1,R2]     
   ])
   #fig.show()

   A = Normal(R1)
   B = Normal(R2)

 
 #  print(A)
 #  print(B)
 #  print(d)

   # Apply separating axis theorem
   cross=[[np.cross(A[i], B[j])/np.linalg.norm(np.cross(A[i], B[j])) for i in range(0,3)] for j in range(0,3)]   

   axis = [A[0], A[1], A[2], B[0], B[1], B[2], cross[0][0], cross[0][1], cross[0][2], cross[1][0], cross[1][1], cross[1][2], cross[2][0], cross[2][1], cross[2][2]]
   
   foundseparate=0

   for ax in axis:
     D= np.abs(np.dot(d, ax))
     D1= np.abs(np.dot(a*A[0], ax))+ np.abs(np.dot(a*A[1], ax))+ np.abs(np.dot(b*A[2], ax))
     D2= np.abs(np.dot(a*B[0], ax))+ np.abs(np.dot(a*B[1], ax))+ np.abs(np.dot(b*B[2], ax))
     foundseparate += D1+D2 < D

  
   noboundaryintersection=0
   for i in range(0,8): 
     noboundaryintersection += (np.sqrt(R1[i][0]*R1[i][0]+R1[i][1]*R1[i][1])<100 or np.sqrt(R2[i][0]*R2[i][0]+R2[i][1]*R2[i][1])<100
         or np.sqrt(R1[i][0]*R1[i][0]+R1[i][1]*R1[i][1])>950 or np.sqrt(R2[i][0]*R2[i][0]+R2[i][1]*R2[i][1])>950
         or R1[i][2]<-2945-200  or R1[i][2]< -2945-200 or R1[i][2]>-2945+200  or R1[i][2]> -2945+200)
         
     

   if np.sum(foundseparate>0 or noboundaryintersection>0):
     return False
   else:
     return True    


noise_level = 0.2

def f(x, noise_level=noise_level, width=a, length=b):
    sum=0
    x=np.array(x).reshape(-1,3)
    print(x)
    for i in x:
        sum+= -1*a**2/np.sum(np.array(i)**2)
    for i in range(0,len(x)):
        for j in range(0,i):
           if(Intersects(np.array(x[i]), np.array(x[j]))):
               sum+= 1e9
    return sum+ np.random.randn() * noise_level


min=[]
max=[]
for i in range(0,110):
  min.append(-950)
  min.append(-950)
  min.append(-2945-200)
  max.append(950)
  max.append(950)
  max.append(-2945+200)


min=np.array(min)
max=np.array(max)

res = gp_minimize(f,                   # the function to minimize
                  [(min[i],max[i])  for i in range(0,330)],      # the bounds on each dimension of x
                  acq_func="EI",       # the acquisition function
                  acq_optimizer="lbfgs", # the acquistion optimizer
                  n_restarts_optimizer=5, # restart optimization with 5 local minima as initial points
                  n_calls=25,          # the number of evaluations of f
                  n_initial_points=25,  # the number of random initialization points
                  noise=0.5**2,        # the noise level (optional)
                  random_state=1234,
                  n_jobs=8)   # the random seed


print(res.x)



from skopt.plots import plot_convergence
plot_convergence(res)






