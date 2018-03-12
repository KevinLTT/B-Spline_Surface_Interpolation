# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import B_Spline_Curve as BCurve
from matplotlib import cm

## generate parameters

def generateParameters( dataPoints, degreeU, degreeV, methodU="uniform", methodV="uniform", alpha=1 ):
    m = dataPoints.shape[0]-1
    n = dataPoints.shape[1]-1
    
    #calculate s0, ..., sm of direction u
    Uij = np.zeros( [m+1,n+1] )
    for j in range( 0, n+1 ):
        Uij[:,j] = BCurve.generateParameters( dataPoints[:,j], degreeU, methodU, alpha )
    s = np.mean( Uij, axis=1 ) 
    
    #calculate t0, ..., tn of direction v
    Vij = np.zeros( [m+1, n+1] )
    for i in range( 0, m+1 ):
        Vij[i,:] = BCurve.generateParameters( dataPoints[i,:], degreeV, methodV, alpha )
    t = np.mean( Vij, axis=0 ) 
    
    return s, t

## generate knots

def generateKnots( dataPoints, degreeU, degreeV, paramU, paramV, methodU="uniform", methodV="uniform" ):
    knotsU = BCurve.generateKnots( dataPoints[:,0], degreeU, paramU, methodU )
    knotsV = BCurve.generateKnots( dataPoints[0], degreeV, paramV, methodV )
    
    return knotsU, knotsV


## calculate control points for each column

def getControlPoints( dataPoints, degreeU, degreeV, paramU, paramV, knotsU, knotsV ):
    m = dataPoints.shape[0] - 1
    n = dataPoints.shape[1] - 1
    
    ctrlPointsColumn = np.zeros( dataPoints.shape )
    for j in range( 0, n+1 ):
        ctrlPointsColumn[:,j] = BCurve.getControlPoints( dataPoints[:,j], degreeU, paramU, knotsU )
    
    ctrlPoints = np.zeros( dataPoints.shape )
    for i in range( 0, m+1 ):
        ctrlPoints[i] = BCurve.getControlPoints( ctrlPointsColumn[i], degreeV, paramV, knotsV )
        
    return ctrlPoints

def BSplineSurfacePoint( ctrlPoints, knotsU, knotsV, degreeU, degreeV, u, v ):
    m = ctrlPoints.shape[0]-1
    n = ctrlPoints.shape[1]-1
    
    point = np.zeros( 3 )
    for i in range( 0, m+1 ):
        for j in range( 0, n+1 ):
            point += ctrlPoints[i][j] * BCurve.BSplineBaseFunction( i, degreeU, m, knotsU, u ) * BCurve.BSplineBaseFunction( j, degreeV, n, knotsV, v)
    
    return point

def BSplineSurface( ctrlPoints, knotsU, knotsV, degreeU, degreeV, piece=100 ):
    x = np.zeros( [piece, piece, 3] )
    for i in range( 0, piece ):
        for j in range( 0, piece ):
            x[i][j] = BSplineSurfacePoint( ctrlPoints, knotsU, knotsV, degreeU, degreeV, i/piece, j/piece )
            
    return x

def BSplineSurfaceInterpolation( dataPoints, degreeU, degreeV, paramMethodU = "uniform", paramMethodV="uniform", knotsMethodU="uniform", knotsMethodV="uniform", piece=100 ):
    paramU, paramV = generateParameters( dataPoints, degreeU, degreeV, paramMethodU, paramMethodV )
    knotsU, knotsV = generateKnots( dataPoints, degreeU, degreeV, paramU, paramV, knotsMethodU, knotsMethodV )
    ctrlPoints = getControlPoints( dataPoints, degreeU, degreeV, paramU, paramV, knotsU, knotsV )
    surface = BSplineSurface( ctrlPoints, knotsU, knotsV, degreeU, degreeV, piece )
    
    return surface, ctrlPoints


