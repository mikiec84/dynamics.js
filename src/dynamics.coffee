# Vector and Matrix extracted from Sylvester.js
class Vector
  constructor: (@els) ->

  # Returns element i of the vector
  e: (i) =>
    return if (i < 1 || i > this.els.length) then null else this.els[i-1]

  # Returns the scalar product of the vector with the argument
  # Both vectors must have equal dimensionality
  dot: (vector) =>
    V = vector.els || vector
    product = 0
    n = this.els.length
    return null if n != V.length
    n += 1
    while --n
      product += this.els[n-1] * V[n-1]
    return product

  # Returns the vector product of the vector with the argument
  # Both vectors must have dimensionality 3
  cross: (vector) =>
    B = vector.els || vector
    return null if this.els.length != 3 || B.length != 3
    A = this.els
    return new Vector([
      (A[1] * B[2]) - (A[2] * B[1]),
      (A[2] * B[0]) - (A[0] * B[2]),
      (A[0] * B[1]) - (A[1] * B[0])
    ])

class Matrix
  constructor: (@els) ->

  # Returns element (i,j) of the matrix
  e: (i,j) =>
    return null if (i < 1 || i > this.els.length || j < 1 || j > this.els[0].length)
    this.els[i-1][j-1]

  # Returns a copy of the matrix
  dup: () =>
    return new Matrix(this.els)

  # Returns the result of multiplying the matrix from the right by the argument.
  # If the argument is a scalar then just multiply all the elements. If the argument is
  # a vector, a vector is returned, which saves you having to remember calling
  # col(1) on the result.
  multiply: (matrix) =>
    returnVector = if matrix.modulus then true else false
    M = matrix.els || matrix
    M = new Matrix(M).els if (typeof(M[0][0]) == 'undefined')
    ni = this.els.length
    ki = ni
    kj = M[0].length
    cols = this.els[0].length
    elements = []
    ni += 1
    while (--ni)
      i = ki - ni
      elements[i] = []
      nj = kj
      nj += 1
      while (--nj)
        j = kj - nj
        sum = 0
        nc = cols
        nc += 1
        while (--nc)
          c = cols - nc
          sum += this.els[i][c] * M[c][j]
        elements[i][j] = sum

    M = new Matrix(elements)
    return if returnVector then M.col(1) else M

  # Returns the transpose of the matrix
  transpose: =>
    rows = this.els.length
    cols = this.els[0].length
    elements = []
    ni = cols
    ni += 1
    while (--ni)
      i = cols - ni
      elements[i] = []
      nj = rows
      nj += 1
      while (--nj)
        j = rows - nj
        elements[i][j] = this.els[j][i]
    return new Matrix(elements)

  # Make the matrix upper (right) triangular by Gaussian elimination.
  # This method only adds multiples of rows to other rows. No rows are
  # scaled up or switched, and the determinant is preserved.
  toRightTriangular: =>
    M = this.dup()
    n = this.els.length
    k = n
    kp = this.els[0].length
    while (--n)
      i = k - n
      if (M.els[i][i] == 0)
        for j in [i + 1...k]
          if (M.els[j][i] != 0)
            els = []
            np = kp
            np += 1
            while (--np)
              p = kp - np
              els.push(M.els[i][p] + M.els[j][p])
            M.els[i] = els
            break
      if (M.els[i][i] != 0)
        for j in [i + 1...k]
          multiplier = M.els[j][i] / M.els[i][i]
          els = []
          np = kp
          np += 1
          while (--np)
            p = kp - np
            # Elements with column numbers up to an including the number
            # of the row that we're subtracting can safely be set straight to
            # zero, since that's the point of this routine and it avoids having
            # to loop over and correct rounding errors later
            els.push(if p <= i then 0 else M.els[j][p] - M.els[i][p] * multiplier)
          M.els[j] = els
    return M

  # Returns the result of attaching the given argument to the right-hand side of the matrix
  augment: (matrix) =>
    M = matrix.els || matrix
    M = new Matrix(M).els if (typeof(M[0][0]) == 'undefined')
    T = this.dup()
    cols = T.els[0].length
    ni = T.els.length
    ki = ni
    kj = M[0].length
    return null if (ni != M.length)
    ni += 1
    while (--ni)
      i = ki - ni
      nj = kj
      nj += 1
      while (--nj)
        j = kj - nj
        T.els[i][cols + j] = M[i][j]

    return T

  # Returns the inverse (if one exists) using Gauss-Jordan
  inverse: =>
    vni = this.els.length
    ki = ni
    M = this.augment(Matrix.I(ni)).toRightTriangular()
    kp = M.els[0].length
    inverse_elements = []
    # Matrix is non-singular so there will be no zeros on the diagonal
    # Cycle through rows from last to first
    ni += 1
    while (--ni)
      i = ni - 1
      # First, normalise diagonal elements to 1
      els = []
      np = kp
      inverse_elements[i] = []
      divisor = M.els[i][i]
      np += 1
      while (--np)
        p = kp - np
        new_element = M.els[i][p] / divisor
        els.push(new_element)
        # Shuffle of the current row of the right hand side into the results
        # array as it will not be modified by later runs through this loop
        if (p >= ki)
          inverse_elements[i].push(new_element)

      M.els[i] = els
      # Then, subtract this row from those above it to
      # give the identity matrix on the left hand side
      for j in [0...i]
        els = []
        np = kp
        np += 1
        while (--np)
          p = kp - np
          els.push(M.els[j][p] - M.els[i][p] * M.els[j][i])
        M.els[j] = els

    return new Matrix(inverse_elements)

  @I = (n) ->
    els = []
    k = n
    n += 1
    while --n
      i = k - n
      els[i] = []
      nj = k
      nj += 1
      while --nj
        j = k - nj
        els[i][j] = if (i == j) then 1 else 0

    new Matrix(els)

# Private Classes
## Dynamics
class Dynamic
  @properties: {}

  constructor: (@options = {}) ->
    for k, v of @options.type.properties
      if !@options[k]? and !v.editable
        @options[k] = v.default

  init: =>
    @t = 0

  next: (step) =>
    @t = 1 if @t > 1
    r = @at(@t)
    @t += step
    r

  at: (t) ->
    [t, t]

class Linear extends Dynamic
  @properties:
    duration: { min: 100, max: 4000, default: 1000 }

  at: (t) ->
    [t, t]

class Gravity extends Dynamic
  @properties:
    bounce: { min: 0, max: 80, default: 40 }
    gravity: { min: 1, max: 4000, default: 1000 }
    expectedDuration: { editable: false }

  constructor: (@options = {}) ->
    @initialForce ?= false
    @options.duration = @duration()
    super @options

  expectedDuration: =>
    @duration()

  duration: =>
    Math.round(1000 * 1000 / @options.gravity * @length())

  bounceValue: =>
    Math.min((@options.bounce / 100), 80)

  gravityValue: =>
    @options.gravity / 100

  length: =>
    bounce = @bounceValue()
    gravity = @gravityValue()
    b = Math.sqrt(2 / gravity)
    curve = { a: -b, b: b, H: 1 }
    if @initialForce
      curve.a = 0
      curve.b = curve.b * 2
    while curve.H > 0.001
      L = curve.b - curve.a
      curve = { a: curve.b, b: curve.b + L * bounce, H: curve.H * bounce * bounce }
    curve.b

  init: =>
    super
    L = @length()
    gravity = @gravityValue() * L * L
    bounce = @bounceValue()

    b = Math.sqrt(2 / gravity)
    @curves = []
    curve = { a: -b, b: b, H: 1 }
    if @initialForce
      curve.a = 0
      curve.b = curve.b * 2
    @curves.push curve
    while curve.b < 1 and curve.H > 0.001
      L = curve.b - curve.a
      curve = { a: curve.b, b: curve.b + L * bounce, H: curve.H * bounce * bounce }
      @curves.push curve

  curve: (a, b, H, t) =>
    L = b - a
    t2 = (2 / L) * (t) - 1 - (a * 2 / L)
    c = t2 * t2 * H - H + 1
    c = 1 - c if @initialForce
    c

  at: (t) =>
    bounce = (@options.bounce / 100)
    gravity = @options.gravity

    i = 0
    curve = @curves[i]
    while(!(t >= curve.a and t <= curve.b))
      i += 1
      curve = @curves[i]
      break unless curve

    if !curve
      v = if @initialForce then 0 else 1
    else
      v = @curve(curve.a, curve.b, curve.H, t)

    [t, v]

class GravityWithForce extends Gravity
  returnsToSelf: true

  constructor: (@options = {}) ->
    @initialForce = true
    super @options

class Spring extends Dynamic
  @properties:
    frequency: { min: 0, max: 100, default: 15 }
    friction: { min: 1, max: 1000, default: 200 }
    anticipationStrength: { min: 0, max: 1000, default: 0 }
    anticipationSize: { min: 0, max: 99, default: 0 }
    duration: { min: 100, max: 4000, default: 1000 }

  at: (t) =>
    frequency = Math.max(1, @options.frequency)
    friction = Math.pow(20, (@options.friction / 100))
    s = @options.anticipationSize / 100
    decal = Math.max(0, s)

    frictionT = (t / (1 - s)) - (s / (1 - s))

    if t < s
      # In case of anticipation
      A = (t) =>
        M = 0.8

        x0 = (s / (1 - s))
        x1 = 0

        b = (x0 - (M * x1)) / (x0 - x1)
        a = (M - b) / x0

        (a * t * @options.anticipationStrength / 100) + b

      yS = (s / (1 - s)) - (s / (1 - s))
      y0 = (0 / (1 - s)) - (s / (1 - s))
      b = Math.acos(1 / A(yS))
      a = (Math.acos(1 / A(y0)) - b) / (frequency * (-s))
    else
      # Normal curve
      A = (t) =>
        Math.pow(friction / 10,-t) * (1 - t)

      b = 0
      a = 1

    At = A(frictionT)

    angle = frequency * (t - s) * a + b
    v = 1 - (At * Math.cos(angle))
    [t, v, At, frictionT, angle]

class SelfSpring extends Dynamic
  @properties:
    frequency: { min: 0, max: 100, default: 15 }
    friction: { min: 1, max: 1000, default: 200 }
    duration: { min: 100, max: 4000, default: 1000 }

  returnsToSelf: true

  at: (t) =>
    frequency = Math.max(1, @options.frequency)
    friction = Math.pow(20, (@options.friction / 100))

    # Normal curve
    A = (t) =>
      1 - Math.pow(friction / 10,-t) * (1 - t)

    At = A(t)
    At2 = A(1-t)

    Ax = (Math.cos(t * 2 * 3.14 - 3.14) / 2) + 0.5
    Ax = Math.pow(Ax, @options.friction / 100)

    angle = frequency * t
    # v = 1 - (At * Math.cos(angle))
    v = Math.cos(angle) * Ax
    [t, v, Ax, -Ax]

class Bezier extends Dynamic
  @properties:
    points: { type: 'points', default: [{x:0,y:0,controlPoints:[{x:0.2,y:0}]},{x:0.5,y:1.2,controlPoints:[{x:0.3,y:1.2},{x:0.8,y:1.2}]},{x:1,y:1,controlPoints:[{x:0.8,y:1}]}] }
    duration: { min: 100, max: 4000, default: 1000 }

  constructor: (@options = {}) ->
    @returnsToSelf = @options.points[@options.points.length - 1].y == 0
    super @options

  B_: (t, p0, p1, p2, p3) ->
    (Math.pow(1 - t, 3) * p0) + (3 * Math.pow(1 - t, 2) * t * p1) + (3 * (1 - t) * Math.pow(t, 2) * p2) + Math.pow(t, 3) * p3

  B: (t, p0, p1, p2, p3) =>
    {
      x: @B_(t, p0.x, p1.x, p2.x, p3.x),
      y: @B_(t, p0.y, p1.y, p2.y, p3.y)
    }

  yForX: (xTarget, Bs) =>
    # Find the right Bezier curve first
    B = null
    for aB in Bs
      if xTarget >= aB(0).x and xTarget <= aB(1).x
        B = aB
      break if B != null

    unless B
      if @returnsToSelf
        return 0
      else
        return 1

    # Find the percent with dichotomy
    xTolerance = 0.0001
    lower = 0
    upper = 1
    percent = (upper + lower) / 2

    x = B(percent).x
    i = 0

    while Math.abs(xTarget - x) > xTolerance and i < 100
      if xTarget > x
        lower = percent
      else
        upper = percent

      percent = (upper + lower) / 2
      x = B(percent).x
      i += 1

    # Returns y at this specific percent
    return B(percent).y

  at: (t) =>
    x = t
    points = @options.points || Bezier.properties.points.default
    Bs = []
    for i of points
      k = parseInt(i)
      break if k >= points.length - 1
      ((pointA, pointB) =>
        B = (t) =>
          @B(t, pointA, pointA.controlPoints[pointA.controlPoints.length - 1], pointB.controlPoints[0], pointB)
        Bs.push(B)
      )(points[k], points[k + 1])
    y = @yForX(x, Bs)
    [x, y]

class EaseInOut extends Dynamic
  @properties:
    friction: { min: 1, max: 1000, default: 500 }
    duration: { min: 100, max: 4000, default: 1000 }

  constructor: (@options = {}) ->
    super
    friction = @options.friction || EaseInOut.properties.friction.default
    points = [
      { x:0, y:0, controlPoints:[{ x:1 - (friction / 1000), y:0 }] },
      { x:1, y:1, controlPoints:[{ x:friction / 1000, y:1 }] }
    ]
    @bezier = new Bezier({
      type: Bezier,
      duration: @options.duration,
      points: points
    })

  at: (t) =>
    @bezier.at(t)

## Helpers
cacheFn = (func) ->
  data = {}
  cachedMethod = ->
    key = ""
    for k in arguments
      key += k.toString() + ","
    result = data[key]
    unless result
      data[key] = result = func.apply(this, arguments)
    result
  cachedMethod

# Browser Support
browserSupportTransform = ->
  browserSupportWithPrefix("transform")

browserSupportPrefixFor = cacheFn (property) ->
  propArray = property.split('-')
  propertyName = ""
  for prop in propArray
    propertyName += prop.substring(0, 1).toUpperCase() + prop.substring(1)
  for prefix in [ "Webkit", "Moz" ]
    k = prefix + propertyName
    if document.body.style[k] != undefined
      return prefix
  ''
browserSupportWithPrefix = (property) ->
  prefix = browserSupportPrefixFor(property)
  return "#{prefix}#{property.substring(0, 1).toUpperCase() + property.substring(1)}" if prefix == 'Moz'
  return "-#{prefix.toLowerCase()}-#{property}" if prefix != ''
  property

# Additional vector tools
lengthVector = (vector) ->
  a = 0
  for e in vector.els
    a += Math.pow(e, 2)
  Math.sqrt(a)
normalizeVector = (vector) ->
  length = lengthVector(vector)
  newElements = []
  for i, e of vector.els
    newElements[i] = e / length
  new Vector(newElements)
combineVector = (a, b, ascl, bscl) ->
  result = []
  for i in [0..2]
    result[i] = (ascl * a.els[i]) + (bscl * b.els[i])
  return new Vector(result)

# Matrix tools
decomposeMatrix = (matrix) ->
  translate = []
  scale = []
  skew = []
  quaternion = []
  perspective = []

  els = matrix.els

  if (els[3][3] == 0)
    return false

  # Normalize the matrix.
  for i in [0..3]
    for j in [0..3]
      els[i][j] /= els[3][3]

  # perspectiveMatrix is used to solve for perspective, but it also provides
  # an easy way to test for singularity of the upper 3x3 component.
  perspectiveMatrix = matrix.dup()

  for i in [0..2]
    perspectiveMatrix.els[i][3] = 0
  perspectiveMatrix.els[3][3] = 1

  # Don't do this anymore, it would return false for scale(0)..
  # if perspectiveMatrix.determinant() == 0
  #   return false

  # First, isolate perspective.
  if els[0][3] != 0 || els[1][3] != 0 || els[2][3] != 0
    # rightHandSide is the right hand side of the equation.
    rightHandSide = new Vector(els[0..3][3])

    # Solve the equation by inverting perspectiveMatrix and multiplying
    # rightHandSide by the inverse.
    inversePerspectiveMatrix = perspectiveMatrix.inverse()
    transposedInversePerspectiveMatrix = inversePerspectiveMatrix.transpose()
    perspective = transposedInversePerspectiveMatrix.multiply(rightHandSide).els

    # Clear the perspective partition
    for i in [0..2]
      els[i][3] = 0
    els[3][3] = 1
  else
    # No perspective.
    perspective = [0,0,0,1]

  # Next take care of translation
  for i in [0..2]
    translate[i] = els[3][i]
    els[3][i] = 0

  # Now get scale and shear. 'row' is a 3 element array of 3 component vectors
  row = []
  for i in [0..2]
    row[i] = new Vector(els[i][0..2])

  # Compute X scale factor and normalize first row.
  scale[0] = lengthVector(row[0])
  row[0] = normalizeVector(row[0])

  # Compute XY shear factor and make 2nd row orthogonal to 1st.
  skew[0] = row[0].dot(row[1])
  row[1] = combineVector(row[1], row[0], 1.0, -skew[0])

  # Now, compute Y scale and normalize 2nd row.
  scale[1] = lengthVector(row[1])
  row[1] = normalizeVector(row[1])
  skew[0] /= scale[1]

  # Compute XZ and YZ shears, orthogonalize 3rd row
  skew[1] = row[0].dot(row[2])
  row[2] = combineVector(row[2], row[0], 1.0, -skew[1])
  skew[2] = row[1].dot(row[2])
  row[2] = combineVector(row[2], row[1], 1.0, -skew[2])

  # Next, get Z scale and normalize 3rd row.
  scale[2] = lengthVector(row[2])
  row[2] = normalizeVector(row[2])
  skew[1] /= scale[2]
  skew[2] /= scale[2]

  # At this point, the matrix (in rows) is orthonormal.
  # Check for a coordinate system flip.  If the determinant
  # is -1, then negate the matrix and the scaling factors.
  pdum3 = row[1].cross(row[2])
  if row[0].dot(pdum3) < 0
    for i in [0..2]
      scale[i] *= -1
      for j in [0..2]
        row[i].els[j] *= -1

  # Get element at row
  rowElement = (index, elementIndex) ->
    row[index].els[elementIndex]

  # Euler angles
  rotate = []
  rotate[1] = Math.asin(-rowElement(0, 2))
  if Math.cos(rotate[1]) != 0
    rotate[0] = Math.atan2(rowElement(1, 2), rowElement(2, 2))
    rotate[2] = Math.atan2(rowElement(0, 1), rowElement(0, 0))
  else
    rotate[0] = Math.atan2(-rowElement(2, 0), rowElement(1, 1))
    rotate[1] = 0

  # Now, get the rotations out
  t = rowElement(0, 0) + rowElement(1, 1) + rowElement(2, 2) + 1.0
  if t > 1e-4
    s = 0.5 / Math.sqrt(t)
    w = 0.25 / s
    x = (rowElement(2, 1) - rowElement(1, 2)) * s
    y = (rowElement(0, 2) - rowElement(2, 0)) * s
    z = (rowElement(1, 0) - rowElement(0, 1)) * s
  else if (rowElement(0, 0) > rowElement(1, 1)) && (rowElement(0, 0) > rowElement(2, 2))
    s = Math.sqrt(1.0 + rowElement(0, 0) - rowElement(1, 1) - rowElement(2, 2)) * 2.0
    x = 0.25 * s
    y = (rowElement(0, 1) + rowElement(1, 0)) / s
    z = (rowElement(0, 2) + rowElement(2, 0)) / s
    w = (rowElement(2, 1) - rowElement(1, 2)) / s
  else if rowElement(1, 1) > rowElement(2, 2)
    s = Math.sqrt(1.0 + rowElement(1, 1) - rowElement(0, 0) - rowElement(2, 2)) * 2.0
    x = (rowElement(0, 1) + rowElement(1, 0)) / s
    y = 0.25 * s
    z = (rowElement(1, 2) + rowElement(2, 1)) / s
    w = (rowElement(0, 2) - rowElement(2, 0)) / s
  else
    s = Math.sqrt(1.0 + rowElement(2, 2) - rowElement(0, 0) - rowElement(1, 1)) * 2.0
    x = (rowElement(0, 2) + rowElement(2, 0)) / s
    y = (rowElement(1, 2) + rowElement(2, 1)) / s
    z = 0.25 * s
    w = (rowElement(1, 0) - rowElement(0, 1)) / s

  quaternion = [x, y, z, w]

  result = {
    translate: translate,
    scale: scale,
    skew: skew,
    quaternion: quaternion,
    perspective: perspective,
    rotate: rotate
  }

  for typeKey, type of result
    for k, v of type
      type[k] = 0 if isNaN(v)

  result

interpolateMatrix = (decomposedA, decomposedB, t) ->
  # New decomposedMatrix
  decomposed = {
    translate: [],
    scale: [],
    skew: [],
    quaternion: [],
    perspective: []
  }

  # Linearly interpolate translate, scale, skew and perspective
  for k in [ 'translate', 'scale', 'skew', 'perspective' ]
    for i in [0..decomposedA[k].length-1]
      decomposed[k][i] = (decomposedB[k][i] - decomposedA[k][i]) * t + decomposedA[k][i]

  # Interpolate quaternion
  qa = decomposedA.quaternion
  qb = decomposedB.quaternion

  angle = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3]

  if angle < 0.0
    for i in [0..3]
      qa[i] = -qa[i]
    angle = -angle

  if angle + 1.0 > .05
    if 1.0 - angle >= .05
      th = Math.acos(angle)
      invth = 1.0 / Math.sin(th)
      scale = Math.sin(th * (1.0 - t)) * invth
      invscale = Math.sin(th * t) * invth
    else
      scale = 1.0 - t
      invscale = t
  else
    qb[0] = -qa[1]
    qb[1] = qa[0]
    qb[2] = -qa[3]
    qb[3] = qa[2]
    scale = Math.sin(piDouble * (.5 - t))
    invscale = Math.sin(piDouble * t)

  for i in [0..3]
    decomposed.quaternion[i] = qa[i] * scale + qb[i] * invscale

  return decomposed

recomposeMatrix = (decomposedMatrix) ->
  matrix = Matrix.I(4)

  # apply perspective
  for i in [0..3]
    matrix.els[i][3] = decomposedMatrix.perspective[i]

  # apply rotation
  quaternion = decomposedMatrix.quaternion
  x = quaternion[0]
  y = quaternion[1]
  z = quaternion[2]
  w = quaternion[3]

  # apply skew
  # temp is a identity 4x4 matrix initially
  skew = decomposedMatrix.skew
  match = [[1,0],[2,0],[2,1]]
  for i in [2..0]
    if skew[i]
      temp = Matrix.I(4)
      temp.els[match[i][0]][match[i][1]] = skew[i]
      matrix = matrix.multiply(temp)

  # Construct a composite rotation matrix from the quaternion values
  matrix = matrix.multiply(new Matrix([[
    1 - 2 * (y * y + z * z),
    2 * (x * y - z * w),
    2 * (x * z + y * w),
    0
  ], [
    2 * (x * y + z * w),
    1 - 2 * (x * x + z * z),
    2 * (y * z - x * w),
    0
  ], [
    2 * (x * z - y * w),
    2 * (y * z + x * w),
    1 - 2 * (x * x + y * y),
    0
  ], [ 0, 0, 0, 1 ]]))

  # apply scale and translation
  for i in [0..2]
    for j in [0..2]
      matrix.els[i][j] *= decomposedMatrix.scale[i]
    matrix.els[3][i] = decomposedMatrix.translate[i]

  matrix

matrixToString = (matrix) ->
  str = 'matrix3d('
  for i in [0..3]
    for j in [0..3]
      str += matrix.els[i][j]
      str += ',' unless i == 3 and j == 3
  str += ')'
  str

transformStringToMatrixString = cacheFn (transform) ->
  matrixEl = document.createElement('div')
  matrixEl.style[browserSupportTransform()] = transform
  document.body.appendChild(matrixEl)
  style = window.getComputedStyle(matrixEl, null)
  result = style.transform || style[browserSupportTransform()]
  document.body.removeChild(matrixEl)
  result

convertToMatrix3d = (transform) ->
  match = transform.match /matrix3?d?\(([-0-9, \.]*)\)/
  if match
    digits = match[1].split(',')
    digits = digits.map(parseFloat)
    if digits.length == 6
      # format: matrix(a, c, b, d, tx, ty)
      elements = [digits[0], digits[1], 0, 0, digits[2], digits[3], 0, 0, 0, 0, 1, 0, digits[4], digits[5], 0, 1]
    else
      elements = digits
  else
    elements = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]

  matrixElements = []
  for i in [0..3]
    matrixElements.push(elements.slice(i * 4,i * 4 + 4))
  new Matrix(matrixElements)

# Private Methods
getFirstFrame = (properties) ->
  frame = {}
  style = window.getComputedStyle(@el, null)
  for k of properties
    v = @el.style[browserSupportWithPrefix(k)]
    v = style[browserSupportWithPrefix(k)] unless v
    frame[k] = v
  frame

parseFrames = (frames) ->
  newFrames = {}
  for percent, properties of frames
    newProperties = {}
    for k, v of properties
      if k != 'transform'
        vString = v + ""
        match = vString.match /([-0-9.]*)(.*)/
        value = parseFloat(match[1])
        unit = match[2]
      else
        value = decomposeMatrix(convertToMatrix3d(transformStringToMatrixString(v)))
        unit = ''
      newProperties[k] = {
        value: value,
        unit: unit
      }
    newFrames[percent] = newProperties
  newFrames

defaultForProperty = (property) ->
  return 1 if property == 'opacity'
  0

animationFrame = (ts) ->
  return if @stopped
  t = 0
  if @ts
    dTs = ts - @ts
    t = dTs / @options.duration
  else
    @ts = ts

  at = @dynamic().at(t)

  animationFrameApply.call(@, at[1], { progress: t })

  if t < 1
    requestAnimationFrame animationFrame.bind(@)
  else
    @animating = false
    @dynamic().init()
    @options.complete?(@)

animationFrameApply = (t, args = {}) ->
  frame0 = @frames[0]
  frame1 = @frames[100]
  progress = args.progress
  progress ?= -1

  transform = ''
  properties = {}
  for k, v of frame1
    value = v.value
    unit = v.unit

    newValue = null
    if progress >= 1
      if @returnsToSelf
        newValue = frame0[k].value
      else
        newValue = frame1[k].value

    if k == 'transform'
      newValue ?= interpolateMatrix(frame0[k].value, frame1[k].value, t)
      matrix = recomposeMatrix(newValue)
      properties['transform'] = matrixToString(matrix)
    else
      unless newValue
        oldValue = null
        oldValue = frame0[k].value if frame0[k]
        oldValue = defaultForProperty(k) unless oldValue?
        dValue = value - oldValue
        newValue = oldValue + (dValue * t)
      properties[k] = newValue

  css(@el, properties)

animationStart = ->
  unless @options.animated
    animationFrameApply.call(@, 1, { progress: 1 })
    return

  @animating = true
  @ts = null
  if @stopped
    @stopped = false
  requestAnimationFrame animationFrame.bind(@)

Animations = []
hasCommonProperties = (props1, props2) ->
  for k, v of props1
    return true if props2[k]?
  false
stopAnimationsForEl = (el, properties) ->
  for animation in Animations
    if animation.el == el and hasCommonProperties(animation.to, properties)
      animation.stop()

# Public Methods
pxProperties = [
  'marginTop', 'marginLeft', 'marginBottom', 'marginRight',
  'paddingTop', 'paddingLeft', 'paddingBottom', 'paddingRight',
  'top', 'left', 'bottom', 'right',
]
css = (el, properties) ->
  for k, v of properties
    unit = ''
    if pxProperties.indexOf(k.toLowerCase()) != -1
      unit = 'px'
    el.style[browserSupportWithPrefix(k)] = "#{v}#{unit}"

# Public Classes
class Animation
  @index: 0

  constructor: (@el, @to, options = {}) ->
    if window['jQuery'] and @el instanceof jQuery
      @el = @el[0]
    @animating = false
    redraw = @el.offsetHeight # Hack to redraw the element
    @frames = parseFrames({
      0: getFirstFrame.call(@, @to),
      100: @to
    })
    @setOptions(options)
    if @options.debugName and Dynamics.InteractivePanel
      Dynamics.InteractivePanel.addAnimation(@)
    Animations.push(@)

  setOptions: (options = {}) =>
    optionsChanged = @options?.optionsChanged

    @options = options
    @options.duration ?= 1000
    @options.complete ?= null
    @options.type ?= Linear
    @options.animated ?= true
    @returnsToSelf = false || @dynamic().returnsToSelf
    @_dynamic = null

    if @options.debugName? and Dynamics.Overrides? and Dynamics.Overrides.for(@options.debugName)
      @options = Dynamics.Overrides.getOverride(@options, @options.debugName)

    @dynamic().init()

    optionsChanged?()

  dynamic: =>
    @_dynamic ?= new @options.type(@options)
    @_dynamic

  start: (options = {}) =>
    options.delay ?= 0
    stopAnimationsForEl(@el, @to)
    if options.delay <= 0
      animationStart.call(@)
    else
      setTimeout animationStart.bind(@), options.delay

  stop: =>
    @animating = false
    @stopped = true

# Export
Dynamics =
  Animation: Animation
  Types:
    Spring: Spring
    SelfSpring: SelfSpring
    Gravity: Gravity
    GravityWithForce: GravityWithForce
    Linear: Linear
    Bezier: Bezier
    EaseInOut: EaseInOut
  css: css

try
  if module
    module.exports = Dynamics
  else
    @Dynamics = Dynamics
catch e
  @Dynamics = Dynamics
