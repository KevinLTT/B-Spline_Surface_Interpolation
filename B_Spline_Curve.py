
# coding: utf-8

# In[50]:


import matplotlib.pyplot as plt
import numpy as np


# In[51]:


## data points and degree

#dataPoints = np.array([
#    [0, 1], [1, 2], [0, 3], [3, 3], [4, 2], [5, 1]#, [6, 2]
#])
#dataPoints = np.array([
#    [0, 1, 1], [1, 2, 0], [0, 3, -1], [3, 3, 0], [4, 2, 0], [5, 1, 0]#, [6, 2]
#])
#degree = 4
#ctrlPoints = np.array([
#    [0, 0], [1, 2], [3, 4], [4, 0]
#])
#plt.plot( dataPoints[:,0], dataPoints[:,1], "b--o" )
#plt.show()


# In[52]:


## parameters

# compute chord length
def chordLength( a, b ):
    dimension = a.shape[0]
    chordLength = 0
    for i in range( 0, dimension ):
        chordLength += (a[i]-b[i])**2
    chordLength = np.sqrt( chordLength)
    return chordLength

# generate parameters
def generateParameters( dataPoints, degree, method="uniform", alpha=1 ):
    n = dataPoints.shape[0]-1
    parameters = np.zeros( n+1 )
    
    if( method == "uniform" ):
        for i in range( 0, n+1 ):
            parameters[i] = i / n
            
    if( method == "chord"):
        Li = np.zeros( n+1 )
        for i in range( 1, n+1 ):
            Li[i] = Li[i-1] + chordLength( dataPoints[i-1], dataPoints[i])
        for i in range( 1, n+1 ):
            parameters[i] = Li[i] / Li[n]
    
    if( method == "centripetal" ):
        Li = np.zeros( n+1 )
        for i in range( 1, n+1 ):
            Li[i] = Li[i-1] + chordLength( dataPoints[i-1], dataPoints[i]) ** alpha
        for i in range( 1, n+1 ):
            parameters[i] = Li[i] / Li[n]
            
    return parameters
          
#parameters = generateParameters( dataPoints, degree, "uniform" )
#print( "parameters: ", parameters)


# In[53]:


## knots

def generateKnots( dataPoints, degree, parameters, method="uniform" ):
    n = dataPoints.shape[0]-1
    m = n + degree + 1
    knots = np.zeros( m+1 )
    for i in range( m-degree, m+1 ):
        knots[i] = 1
    
    if( method == "uniform" ):
        for j in range( 1, m-2*degree ):
            knots[j+degree] = j / ( n - degree + 1)
        
    if( method == "average"):
        for j in range( 1, m-2*degree ):
            tSum = 0
            for i in range( j, j+degree ):
                tSum += parameters[i]
            knots[j+degree] = tSum / degree
    
    return knots

#knots = generateKnots( dataPoints, degree, parameters, "uniform" )
#print( "knots: ", knots )


# In[54]:


## B-Spline Base Function

# n+1: number of control points 
# k: degree = order - 1
# knots = { ti: i = 0, 1, ..., n+k+1 }
def BSplineBaseFunction( i, k, n, knots, t ):
    if( k == 0 ):
        if t >= knots[i] and t < knots[i+1] :
            return 1
        else:
            return 0
    
    factor1 = 0
    if i+k <= n+k+1 and knots[i+k]-knots[i] != 0:
        factor1 = ( t - knots[i] ) / ( knots[i+k] - knots[i] )
    term1 = factor1 * BSplineBaseFunction( i, k-1, n, knots, t )
    
    factor2 = 0
    if i+k+1 <= n+k+1 and knots[i+k+1]-knots[i+1] != 0:
        factor2 = ( knots[i+k+1] - t ) / ( knots[i+k+1] - knots[i+1] )
    term2 = 0
    if i+1 <= n+k+1:
        term2 = factor2 * BSplineBaseFunction( i+1, k-1, n, knots, t )
    
    return term1 + term2 

def drawBaseFunction( degree, n, knots, piece=100 ):
    x = np.arange( 0, 1, 1 / piece )
    Base = np.zeros( [n+1, piece] )
    for i in range( 0, n+1 ):    
        for j in range( 0, piece ):
            Base[i][j] = BSplineBaseFunction( i, degree, n, knots, j/piece )
        plt.plot( x, Base[i])
    plt.show()
    
#n = dataPoints.shape[0]-1
#drawBaseFunction( degree, n, knots )

#a = BSplineBaseFunction( 5, 3, n, knots, 1)
#print( a )


# In[55]:


## compute control points

def getControlPoints( dataPoints, degree, parameters, knots ):
    n = dataPoints.shape[0] - 1
    dimension = dataPoints.shape[1]
    
    N = np.zeros( [n+1, n+1] )
    for i in range( 0, n+1 ):
        for j in range( 0, n+1 ):
            N[i][j] = BSplineBaseFunction( j, degree, n, knots, parameters[i] )
    
    N[n][n] = 1
    Ninverse = np.linalg.inv(N)
    
    P = Ninverse.dot( dataPoints )
    return P

#ctrlPoints = getControlPoints( dataPoints, degree, parameters, knots )
#plt.plot( ctrlPoints[:,0], ctrlPoints[:,1], "r--o" )
#plt.show()


# In[56]:


def BSplinePoint( ctrlPoints, knots, degree, t ):
    n = ctrlPoints.shape[0] - 1
    dimension = ctrlPoints.shape[1]
    point = np.zeros( dimension )
    for i in range( 0, n+1 ):
        for j in range ( 0, dimension ):
            point[j] += ctrlPoints[i][j] * BSplineBaseFunction( i, degree, n, knots, t )
            
    return point

def BSplineCurve( ctrlPoints, knots, degree, piece=100 ):
    dimension = ctrlPoints.shape[1]
    x = np.zeros( [piece, dimension] )
    for i in range( 0, piece ):
        x[i] = BSplinePoint( ctrlPoints, knots, degree, i/piece )
    
    return x


#piece = 100
#curve = BSplineCurve( ctrlPoints, knots, degree, piece )
#plt.plot( curve[:,0], curve[:,1] )
#plt.plot( dataPoints[:,0], dataPoints[:,1], "r--o" )
#plt.plot( ctrlPoints[:,0], ctrlPoints[:,1], "y--o" )
#plt.show()


# In[57]:


def BSplineCurveInterpolation( dataPoints, degree, parameterMethod, knotsMethod, piece=100 ):
    parameters = generateParameters( dataPoints, degree, parameterMethod )
    knots = generateKnots( dataPoints, degree, parameters, knotsMethod )
    ctrlPoints = getControlPoints( dataPoints, degree, parameters, knots )
    curve = BSplineCurve( ctrlPoints, knots, degree, piece )
    
    return curve, ctrlPoints


#dataPoints = np.array( [
#    [ [0, 0, 60], [0, 0, 30], [0, 0, 0], [3, 20, 0], [6, 40, -5] ],
#    [ [20, 0, 60], [20, 7, 30], [20, 7, 0], [23, 20, 0], [26, 43, -5] ],
#    [ [40, 0, 60], [40, 7, 30], [40, 7, 0], [43, 20, 0], [46, 43, -5] ],
#    [ [60, 0, 60], [60, 0, 30], [60, 0, 0], [63, 20, 0], [66, 40, -5] ]
#])

#curve, ctrlPoints = BSplineCurveInterpolation( dataPoints, degree, parameterMethod="uniform", knotsMethod="uniform" )
#print( dataPoints[0] )
#curve, ctrlPoints = BSplineCurveInterpolation( dataPoints[0], 3, parameterMethod="uniform", knotsMethod="uniform" )
#plt.plot( dataPoints[0,:,0], dataPoints[0,:,1], "r--o" )
#plt.plot( ctrlPoints[:,0], ctrlPoints[:,1], "y--o" )
#plt.plot( curve[:,0], curve[:,1] )
#plt.show()


