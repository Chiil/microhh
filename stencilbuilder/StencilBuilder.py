#!/usr/bin/python

import numpy as np

# Define the location arrays
uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
sloc = np.array([0,0,0])
zloc  = 0
zhloc = 1

# Check locations
def checkLocs(left, right):
  loc = np.array([ None, None, None ])
  for i in range(len(loc)):
    if (left.loc[i] == None and right.loc[i] == None):
      continue
    elif (left.loc[i] == None):
      loc[i] = right.loc[i]
    elif (right.loc[i] == None):
      loc[i] = left.loc[i]
    elif (left.loc[i] != right.loc[i]):
      raise (Exception("Left and right are not at the same grid location"))
    else:
      loc[i] = left.loc[i]

  return loc

# Base Node class
class Node(object):
  def __add__(self, right):
    return NodeOperator(self, right, "+")

  def __sub__(self, right):
    return NodeOperator(self, right, "-")

  def __mul__(self, right):
    return NodeOperator(self, right, "*")

class NodeOperator(Node):
  def __init__(self, left, right, operatorString):
    self.left = left
    self.right = right
    self.operatorString = operatorString
    self.depth = max(left.depth, right.depth)
    if (self.depth > 1):
      self.pad = 2
    else:
      self.pad = 0

    self.loc = checkLocs(left, right)

  def getString(self, i, j, k, pad, plane):
    ob = '( '
    cb = ' )'

    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(1, self.depth):
        lb = lb + '\n'

      return "{ob}{0}\n{lb}{ws}{os} {1}{cb}".format(self.left.getString(i, j, k, pad, plane),
                                                    self.right.getString(i, j, k, pad, plane),
                                                    ws=ws, lb=lb, ob=ob,
                                                    os=self.operatorString, cb=cb)
    else:
      return "{ob}{0} {os} {1}{cb}".format(self.left.getString(i, j, k, pad, plane),
                                           self.right.getString(i, j, k, pad, plane),
                                           ob=ob, os=self.operatorString, cb=cb)

class NodeStencilFour(Node):
  def __init__(self, inner, dim, c0, c1, c2, c3):
    self.inner = inner
    self.depth = inner.depth + 1
    self.pad = 6

    self.dim = dim
    self.loc = np.copy(inner.loc)
    self.loc[dim] = not self.loc[dim]

    self.c0 = c0
    self.c1 = c1
    self.c2 = c2
    self.c3 = c3

  def getString(self, i, j, k, pad, plane):
    i0 = i1 = i2 = i3 = i
    j0 = j1 = j2 = j3 = j
    k0 = k1 = k2 = k3 = k
    if (self.dim == 0):
      i0 += -1-self.loc[0]
      i1 +=   -self.loc[0]
      i2 += +1-self.loc[0]
      i3 += +2-self.loc[0]
    elif (self.dim == 1):
      j0 += -1-self.loc[1]
      j1 +=   -self.loc[1]
      j2 += +1-self.loc[1]
      j3 += +2-self.loc[1]
    elif (self.dim == 2):
      k0 += -1-self.loc[2]
      k1 +=   -self.loc[2]
      k2 += +1-self.loc[2]
      k3 += +2-self.loc[2]

    ob = '( '
    cb = ' )'

    newplane = np.copy(plane)
    if (self.depth-1 < 2):
      newplane[self.dim] = 1

    if (self.depth > 1):

      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "{ob}{c0}*{0}\n{lb}{ws}+ {c1}*{1}\n{lb}{ws}+ {c2}*{2}\n{lb}{ws}+ {c3}*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad, newplane),
          self.inner.getString(i1, j1, k1, pad, newplane),
          self.inner.getString(i2, j2, k2, pad, newplane),
          self.inner.getString(i3, j3, k3, pad, newplane),
          ws=ws, lb=lb, ob=ob, cb=cb,
          c0=self.c0, c1=self.c1, c2=self.c2, c3=self.c3)
    else:
      return "{ob}{c0}*{0} + {c1}*{1} + {c2}*{2} + {c3}*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad, newplane),
          self.inner.getString(i1, j1, k1, pad, newplane),
          self.inner.getString(i2, j2, k2, pad, newplane),
          self.inner.getString(i3, j3, k3, pad, newplane),
          ob=ob, cb=cb,
          c0=self.c0, c1=self.c1, c2=self.c2, c3=self.c3)

def formatIndex(n, nstr):
  if (n > 0):
    nn = "+{0}{1}".format(nstr, abs(n))
  elif (n < 0):
    nn = "-{0}{1}".format(nstr, abs(n))
  else:
    nn = " " * (len(nstr) + 2)
  return nn

class Field(Node):
  def __init__(self, name, loc):
    self.name  = name
    self.depth = 0
    self.loc   = np.copy(loc)

  def getString(self, i, j, k, pad, plane):
    compact = True
    ii = jj = kk = ""
    if (plane[0] or not compact):
      ii = formatIndex(i, "ii")
    if (plane[1] or not compact):
      jj = formatIndex(j, "jj")
    if (plane[2] or not compact):
      kk = formatIndex(k, "kk")
    return "{0}[ijk{1}{2}{3}]".format(self.name, ii, jj, kk)

# Vector class representing a profile
class Vector(Node):
  def __init__(self, name, loc):
    self.name  = name
    self.depth = 0
    self.loc   = np.array([ None, None, loc ])

  def getString(self, i, j, k, pad, plane):
    kk = formatIndex(k, "")
    return "{0}[k{1}]".format(self.name, kk)

class Scalar(Node):
  def __init__(self, name):
    self.name  = name
    self.depth = 0
    self.loc   = np.array([ None, None, None ])

  def getString(self, i, j, k, pad, plane):
    return "{0}".format(self.name)

# Define functions.
def interpx(inner):
  return NodeStencilFour(inner, 0, "ci0", "ci1", "ci2", "ci3")
def interpy(inner):
  return NodeStencilFour(inner, 1, "ci0", "ci1", "ci2", "ci3")
def interpz(inner):
  return NodeStencilFour(inner, 2, "ci0", "ci1", "ci2", "ci3")

def gradx(inner):
  return NodeStencilFour(inner, 0, "cg0", "cg1", "cg2", "cg3")
def grady(inner):
  return NodeStencilFour(inner, 1, "cg0", "cg1", "cg2", "cg3")
def gradz(inner):
  return NodeStencilFour(inner, 2, "cg0", "cg1", "cg2", "cg3")

def printStencil(lhs, rhs, operator):
  checkLocs(lhs, rhs)

  index = "[ijk]"
  indent = len(lhs.name) + len(index) + len(operator) + 2

  plane = np.array([0,0,0])
  print("{0}{1} {2} {3};".format(lhs.name, index, operator, rhs.getString(0, 0, 0, indent, plane)))

