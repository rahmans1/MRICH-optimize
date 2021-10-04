from skopt import gp_minimize
from skopt.plots import plot_convergence

def Normal(R):
   A1 = np.cross(R[7]-R[2], R[3]-R[2])
   A1 /= np.linalg.norm(A1)
   A2 = np.cross(R[6]-R[1], R[2]-R[1])
   A2 /= np.linalg.norm(A2)
   A3 = np.cross(R[5]-R[6], R[7]-R[6])
   A3 /= np.linalg.norm(A3)

   print(R)
   print(A1)
   print(A2)
   print(A3)

   print(np.dot(A1,A3))
   print(np.dot(A1,A2))
   print(np.dot(A2,A3))
   
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
   fig.show()

   A = Normal(R1)
   B = Normal(R2)

 
   print(A)
   print(B)
   print(d)

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




     
