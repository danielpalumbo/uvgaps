"""
Routines for finding circular gaps in u-v coverage.
Written by Daniel Palumbo and Linus Hamilton.
"""

import ehtim as eh
import matplotlib.pyplot as plt
import numpy
import scipy.spatial
import numpy as np

def largestEmptyCircle(xs, ys):
    """ Finds the largest empty circle, among the points (x,y),
    which cannot freely slide without hitting a point. """
    points = numpy.stack([xs, ys], axis=1)
    triangulation = scipy.spatial.Delaunay(points)
    triangleIndices = triangulation.simplices
    # `triangleIndices` is an array with shape (?,3)
    # whose rows are (i,j,k) indices of the vertices
    # (points[i], points[j], points[k]) in each Delaunay triangle.
    triangles = points[triangleIndices] # Wow, numpy indexing is cool
    # Remove strictly obtuse triangles, since they can always "slide"
    nonobtuseTriangles = removeObtuse(triangles)
    # Each Delaunay triangle begets a circumcircle containing no other points.
    # Now we just find the biggest one.
    circleRadii = circumradii(nonobtuseTriangles)

    biggestCircleIndex = numpy.argmax(circleRadii)
    return (circumcenter(nonobtuseTriangles[biggestCircleIndex]),
            circleRadii[biggestCircleIndex])

def getTriangles(xs, ys):
    """ Finds the largest empty circle, among the points (x,y),
    which cannot freely slide without hitting a point. """
    points = numpy.stack([xs, ys], axis=1)
    triangulation = scipy.spatial.Delaunay(points)
    triangleIndices = triangulation.simplices
    # `triangleIndices` is an array with shape (?,3)
    # whose rows are (i,j,k) indices of the vertices
    # (points[i], points[j], points[k]) in each Delaunay triangle.
    triangles = points[triangleIndices] # Wow, numpy indexing is cool
    # Remove strictly obtuse triangles, since they can always "slide"
    nonobtuseTriangles = removeObtuse(triangles)
    return nonobtuseTriangles

def circumcenter(triangle):
    """ `triangle`: an array of shape (3,2)
    returns an array of shape (2,): the incenter of the triangle. """
    # get side lengths
    A,B,C = triangle
    a,b,c = map(numpy.linalg.norm, (B-C, C-A, A-B))
    # use barycentric coordinate formula
    aw,bw,cw = (a**2*(b**2+c**2-a**2),
                b**2*(c**2+a**2-b**2),
                c**2*(a**2+b**2-c**2))
    return (aw*A + bw*B + cw*C)/(aw+bw+cw)

def removeObtuse(triangles):
    """ `triangles`: an array of shape (?,3,2)
    returns an array of shape (something,3,2) removing obtuse triangles. """
    # get side lengths
    A,B,C = triangles[:,0,:], triangles[:,1,:], triangles[:,2,:]
    a,b,c = [numpy.linalg.norm(vecs, axis=1) for vecs in (A-B,B-C,C-A)]
    # Give a little buffer, so we don't accidentally remove right triangles
    epsilon = 1e-6
    # A triangle is obtuse if (biggest side)^2 is bigger than the sum of
    # the squares of the other two sides
    notObtuse = (a**2+b**2-c**2 > -epsilon * (a+b+c)**2) &\
                (b**2+c**2-a**2 > -epsilon * (a+b+c)**2) &\
                (c**2+a**2-b**2 > -epsilon * (a+b+c)**2)
    return triangles[notObtuse]
    
def circumradii(triangles):
    """ `triangles`: an array of shape (?,3,2)
    returns an array of shape (?,) with the circumradii of the triangles. """
    # get side lengths
    A,B,C = triangles[:,0,:], triangles[:,1,:], triangles[:,2,:]
    a,b,c = [numpy.linalg.norm(vecs, axis=1) for vecs in (A-B,B-C,C-A)]
    # use formula
    return a*b*c/((a+b+c)*(a+b-c)*(b+c-a)*(c+a-b))**0.5

def find_gaps(u, v):
    """
    Given a set of u,v points, return the centers and radii of all gaps
    """

    triangles = getTriangles(u, v)
    triangle_radii = np.sqrt(np.sum(triangles**2,axis=2))
    triangle_radii[triangle_radii==0] = 1e100

    circleRadii = circumradii(triangles)
    coords = [circumcenter(triangle) for triangle in triangles]
    gap_centers = np.array(coords)
    radii = [np.sqrt(coord[0]**2+coord[1]**2) for coord in coords]
    gap_radii = np.array(circleRadii)
    return gap_centers, gap_radii
