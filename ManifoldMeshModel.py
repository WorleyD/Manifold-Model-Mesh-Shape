# SETUP + CONFIG 
import numpy as np
import geomproc 
import os
import sys
import math
from math import cos, sin, atan, acos
import time
import scipy.linalg

PREC = 5

def roundMatrix(M, prec=PREC):
	# round a matrix's elements to the given precision
	for i in range(M.shape[0]):
		for j in range(M.shape[1]):
			M[i][j] = round(M[i][j], prec)

	return M

def roundVector(V, prec=PREC):
	for i in range(len(V)):
		V[i] = round(V[i], prec)
	return V

def matrixEq(m1, m2):
	if m1.shape != m2.shape:
		return False

	for i in range(m1.shape[0]):
		for j in range(m1.shape[1]):
			if round(m1[i][j], 5) != round(m2[i][j], PREC):
				return False
	return True

# translate the triangle at index i to the origin
# Assumes the triangle has been canonicalized so the translation exists in the mesh
def applyTransformation(mesh, index):
	v[0] = mesh.translation[i]
	triangle = np.array([list(mesh.vertex[i]) for i in face])

	for i in range(3):
			for j in range(3):
				triangle[i][j] -= v0[i]

	return triangle

def getTrace(m):
	s = 0
	if m.shape[0] != m.shape[1]:
		print("Cannot take trace of non-square matrix \n", m)
		sys.exit()
	for i in range(m.shape[0]):
		s += round(m[i][i], PREC)
	return s


#matrix exp
def m_exp(m):
	u = m[0][1]
	v = m[1][1]
	if v == 0:
		return np.array([[1,u],[0,1]])
	
	return np.array([[1, (u/v)*(math.exp(v)-1)], [0, math.exp(v)]])


#matrix log
def m_log(m):
	U = m[0][1]
	V = m[1][1]
	v = math.log(abs(V))
	if V == 1:
		u = U
	else:
		u = U*v/(V-1)
	return np.array([[1,u],[0,v]])

# matrix log for 3x3 rotation matrices
def rot_log(m):
	trace = getTrace(m)
	if trace < -1:
		#print("adjusting trace of a rotation matrix due to floating point errors")
		#print("trace before:", trace)
		trace -= (trace + 1)
		#print("trace after:", trace)
	if trace > 3:
		#print("adjusting trace of a rotation matrix due to floating point errors")
		#print("trace before:", trace)
		trace -= (trace-3)
		#print("trace after:", trace)
	
	theta = acos((trace-1)/2)
	if theta == 0:
		return [[0 for _ in range(m.shape[0])] for _ in range(m.shape[1])]
	else:
		return (theta/(2*sin(theta)))*(m - np.matrix.transpose(m))
#sample the geodesic path from p to q at t
# t is between 0 and 1 inclusive
def geodesic_GA(p,q,t):
	return np.matmul(p, m_exp(t* m_log(np.matmul(np.linalg.inv(p), q))))

def geodesic_Rot(p,q,t):
	a = np.matmul(np.linalg.inv(p), q)
	b = t* rot_log(a)
	if t == 0:
		b = np.array([[1,0,0],[0,1,0],[0,0,1]])
	c =  scipy.linalg.expm(b)
	d = np.matmul(p, c)
	return d

def geodesic_Trans(p,q,t):
	return np.array([geodesic_Scale(p[i],q[i], t) for i in range(3)])

def geodesic_Scale(p,q,t):
	return (1-t)*p + t*q

def frob_squared(m):
	#frobenius norm of M
	t = 0
	for row in m:
		for e in row:
			t+= e*e
	return t

#takes two (R,A,S) tuples and calculates the distance between them squared
def d_Gt_squared(t1, t2):
	dist = frob_squared(rot_log(np.matmul(np.matrix.transpose(t1[0]), t2[0])))
	if np.isnan(dist):
		print("Issue calculating distance")
		print("tuple 1:\n", t1[0],"\n",t1[1],"\n",t1[2])
		print("tuple 2:\n", t2[0],"\n",t2[1],"\n",t2[2])
		sys.exit()
	dist += frob_squared(m_log(np.matmul(np.linalg.inv(t1[1]), t2[1])))
	dist += abs(math.log(t2[2]/t1[2]))**2
	return dist

def d_manifold_squared(m1, m2):
	#distance between two points on the manifold
	t = 0
	for i in range(len(m1.canonical)): # This is len of N, as are all arrays listed below 
		gt1 = (m1.rot[i], m1.ga[i], m1.scale[i])
		gt2 = (m2.rot[i], m2.ga[i], m2.scale[i])
		t += d_Gt_squared(gt1,gt2)
	return t

# can1 = [p1x, p2x] can2 = [p1y, p2y]
# p2x = [x2x, y2x] p2y = [x2y,y2y]
def getGa(can1, can2):
	
	p2x = [can1[i][1] for i in range(2)]
	p2y = [can2[i][1] for i in range(2)]

	V = p2y[1]/(p2x[1])
	U = (p2y[0] - p2x[0])/(p2x[1])
	n = np.array([[1, U],[0, V]])
	return n

#scaling factor for the triangle 
def getScale(t1, t2):
	return np.linalg.norm(t2)/np.linalg.norm(t1)


def rotation_matrix_from_vectors_3d(a, b):
    b = b / np.linalg.norm(b) # normalize a
    a = a / np.linalg.norm(a) # normalize b
    v = np.cross(a, b)
    # s = np.linalg.norm(v)
    c = np.dot(a, b)

    v1, v2, v3 = v
    h = 1 / (1 + c)

    Vmat = np.array([[0, -v3, v2],
                  [v3, 0, -v1],
                  [-v2, v1, 0]])

    R = np.eye(3, dtype=np.float64) + Vmat + (Vmat.dot(Vmat) * h)
    return R

def rotation_matrix_from_vectors_2d(a, b):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    b = b / np.linalg.norm(b) # normalize a
    a = a / np.linalg.norm(a) # normalize b

    mar =  [[a[0]*b[0] + a[1]*b[1],   b[0]*a[1]-a[0]*b[1]], 
    		[a[0]*b[1] - b[0]*a[1],   a[0]*b[0] + a[1]*b[1]]]
    return np.array(mar)

def canonicalize(mesh):
	for face in mesh.face:
		# Get 3x3 triangle matrix
		triangle = np.array([list(mesh.vertex[i]) for i in face])

		# WLOG we can translate triangle so v0 is at the origin
		# So get that translation
		v0 = [triangle[i][0] for i in range(3)]
		#Store the translation just in case
		mesh.translation.append(v0)

		# Apply transformation to the triangle
		for i in range(3):
			for j in range(3):
				triangle[i][j] -= v0[i]

		#isolate columns in edge matrix for triangle
		p1 = np.array([triangle[i][1] for i in range(3)])
		p2 = np.array([triangle[i][2] for i in range(3)])

		#solve for rotation matrix (from p1 to [||p1||, 0, 0])
		R = rotation_matrix_from_vectors_3d(p1, np.array([np.linalg.norm(p1), 0, 0]))
		# R is a matrix that takes a triangle with v1=[0,0,0], v2=[x,y,z], v3=[something] and makes v2=[norm(v2),0,0]
		# This also rotates v3, from here. We need to keep v2 fixed and rotate v3 so that v3.z = 0
		
		Rt = np.matmul(R, triangle)
		#ensure we actually have zeros where we want to just in case of floating point error
		Rt = roundMatrix(Rt, 10)

		y = Rt[1][2]
		z = Rt[2][2]
		v2 = np.array(y,z)
		#find c1 and c2 
		R2 = rotation_matrix_from_vectors_2d(np.array([y, z]), np.array([math.sqrt(y**2 + z**2),0]))

		c1 = R2[0][0]
		c2 = R2[1][0]
		n = np.array([[1,0,  0],
					  [0, c1, -1*c2],
					  [0, c2, c1]])

		Rx = np.matmul(n, R)	
		
		R_fin = np.matmul(Rx, R)
		
		res = np.matmul(Rx, triangle)
		if round(res[1][1], 5) != 0 or round(res[2][1], 5) != 0 or round(res[2][2]) != 0:
			print("Canonicalization failed. This shouldn't happen")
			print("\'canonical\' triangle was: ")
			print(res)
			print("From originating triangle:\n", triangle)
			print("With Rx:\n", Rx)
			print("And R:\n", R)
			print("Exiting...")
			sys.exit()
		#add canonical triangle to the mesh
		R_fin_2d = np.array([[res[0][1], res[0][2]],[res[1][1], res[1][2]]])
		mesh.rot.append(roundMatrix(R_fin))
		mesh.canonical.append(R_fin_2d)


def redefine(mesh, base):
	# model mesh as a deformation of base
	# ASSUMPTIONS: There is a one to one correspondence between triangles, and they are indexed in correspondence 

	for index in range(N):
		if index >= len(mesh.canonical):
			break

		mesh_triangle = mesh.canonical[index]
		base_triangle = base.canonical[index]
		#After canonicalizing both triangles are of the form [x1, x1]
		#													 [0,  y2]
		S = base_triangle[0][0]/mesh_triangle[0][0]
		mesh.scale.append(S)
		A = getGa(S*mesh_triangle, base_triangle)
		mesh.ga.append(A)

		res = np.matmul(mesh.ga[index], S*mesh_triangle)
		if not matrixEq(res, base_triangle):
			print("Ga not correct")
			print("Canonical template triangle is:\n", base_triangle)
			print("Canonical mesh triangle (t) is:\n", mesh_triangle)
			print("Scale between triangles:", S)
			print("GA:\n", mesh.ga[index])
			print("GA*S*t:\n", res)
			sys.exit()
		
	#When done each triangle will have its R A S tuple stores in lists of the mesh

def getTriangleFromCan(can, rot, trans, scale):
	# Dont think we need scale but keeping it in case
	#can = [x1, x2]
	#      [0,  y2] so put it in 3D matrix
	x1 = can[0][0]
	x2 = can[0,1]
	y2 = can[1][1] 
	#form we want: [0,x1,x2]
	#              [0, 0,y2]
	#              [0, 0, 0]

	newcan = np.array([[0,x1,x2], [0, 0, y2], [0,0,0]])
	t1 = np.matmul(np.linalg.inv(rot), newcan)
	for i in range(3):
		for j in range(3):
			t1[i][j] += trans[i]
	return t1 



if __name__ == "__main__":

	RUN_TESTS = True
	LOAD_MODEL = True
	TEST_PATHS = False

	#### METHOD TESTS
	if RUN_TESTS:
		print("Running tests:")
		print()
		print("Testing Matrix log and exp")
		M = np.array([[1, .7],[0, .33]])
		print("M:\n",M)
		print()
		print("exp(M):\n",m_exp(M))
		print()
		print("log(M):\n",m_log(M))
		print()
		print("exp(log(M)):\n",m_exp(m_log(M)))
		print()
		print("log(exp(M)):\n", m_log(m_exp(M)))

		print()
		print("Testing geodesic path")
		t1 = np.array([[1, 0.5], [0, 3]])
		t2 = np.array([[1, 4], [0, 7]])
		print("p1:\n",t1)
		print("p2:\n", t2)
		t = -0.25
		while t <= 1.4:
			print("path at t={}".format(t))
			print(geodesic_GA(t1,t2,t))
			t += 0.25
			print()

		print("Testing mapping canonical triangles")
		t1 = np.array([[7,3],[0,4],[0,0]])
		t2 = np.array([[1,-3],[0,4], [0,0]])
		print("t1:\n", t1)
		print("t2:\n", t2)
		S = t2[0][0]/t1[0][0]
		t1 = S*t1
		print("Scale between triangles:", S)

		g = getGa(t1,t2)
		print("G_A:\n", g)
		t12d = np.array([[7,3],[0,4]])
		print("t1*G_A=t2:\n", S*np.matmul(g, t12d))

	#### CONSTRUCT MODEL
	if LOAD_MODEL:
		start = time.time()
		DIR = "toscahires-asci/cat/" 
		FILE = DIR + ""

		TEMPLATE = "cat0.tri"

		template = geomproc.load(DIR+TEMPLATE)
		template.filename = TEMPLATE
		N = len(template.face) #number of triangles
		print("Template with {} triangles loaded. Canonicalizing".format(N))
		canonicalize(template)
		I2 = np.array([[1,0],[0,1]])
		template.ga = [I2 for _ in range(N)]
		template.scale = [1 for _ in range(N)]

		print("Template Canonicalized. Remeshing samples")

		meshes = [template]
		count = 1
		for filename in os.listdir(DIR):
			# Need to switch to .obj to use geomproc library
			if filename == TEMPLATE:
				continue # dont want to reprocess the template
			
			if not filename.endswith(".obj") and not filename.endswith(".tri"):
				continue

			print("loading, canonicalizing, and redefining",filename)
			mesh = geomproc.load(DIR + filename)
			mesh.filename = filename
			if N != len(mesh.face):
				print("Size mismatch in meshes. Template has {} faces and current mesh has {}. Exiting...".format(N, len(mesh.face)))
				sys.exit()

			# Finds the R matrices to canonicalize each triangle
			canonicalize(mesh) 
			# Finds the A and S matrices to scale and convert each
			# canonical triangle to the corresponding base triangle
			redefine(mesh, template) 
			meshes.append(mesh)
			count+=1

		print("Samples remeshed")
		end = time.time()
		print("It took {} seconds to load a model with {} meshes, each with {} triangles".format(end-start, count, N))


		#### RESULTS
		m1 = meshes[1]
		m2 = meshes[2]
		m3 = meshes[3]
		print("Geodesic Distance Tests:")
		print("The distance between mesh {} and mesh {} is: {}".format(m1.filename, m2.filename, math.sqrt(d_manifold_squared(m1, m2))))
		print("The distance between mesh {} and mesh {} is: {}".format(m1.filename, m3.filename, math.sqrt(d_manifold_squared(m1, m3))))
		print("The distance between mesh {} and mesh {} is: {}".format(m2.filename, m3.filename, math.sqrt(d_manifold_squared(m2, m3))))
		print("The distance between mesh {} and mesh {} is: {}".format(m3.filename, m3.filename, math.sqrt(d_manifold_squared(m3, m3))))


		if TEST_PATHS:
			OUTPUT_DIR=DIR + "out/"
			print("Geodesic Path Test:")
			m1 = meshes[1]
			m2 = meshes[4]
			verts1 = m1.vertex
			verts2 = m2.vertex
			print("Sampling the geodesic path between {} and {}".format(m1.filename, m2.filename))

			print("Interpolation:")
			t = 0
			while t <= 1:
				mesh = geomproc.mesh()
				mesh.vertex = verts1
				newTriangles = []
				print("Calculating mesh at t = {}".format(t))
				for i in range(N):
					c1 = m1.canonical[i]
					c2 = m2.canonical[i]
					r1 = m1.rot[i]
					r2 = m2.rot[i]
					t1 = m1.translation[i]
					t2 = m2.translation[i]
					s1 = m1.scale[i]
					s2 = m2.scale[i] 
					inter_c = geodesic_GA(c1,c2,t)
					inter_r = geodesic_Rot(r1,r2,t)
					inter_t = geodesic_Trans(t1,t2,t)
					inter_s = geodesic_Scale(s1,s2,t)
					newTriangles.append(getTriangleFromCan(inter_c, inter_r, inter_t, inter_s))

				mesh.face = np.array(newTriangles)
				mesh.save(OUTPUT_DIR+"path_at_{}.obj".format(t))	
				print("mesh \'path_at_{}.obj\' saved.".format(t))
				t+= 0.25

			print("Extrapolation:")
			print("Calculating mesh at t=-0.25")
			t = -1*0.25
			mesh = geomproc.mesh()
			mesh.vertex = verts1
			newTriangles = []
			for i in range(N):
				c1 = m1.canonical[i]
				c2 = m2.canonical[i]
				r1 = m1.rot[i]
				r2 = m2.rot[i]
				t1 = m1.translation[i]
				t2 = m2.translation[i] 
				s1 = m1.scale[i]
				s2 = m2.scale[i]
				inter_c = geodesic_GA(c1,c2,t)
				inter_r = geodesic_Rot(r1,r2,t)
				inter_t = geodesic_Trans(t1,t2,t)
				inter_s = geodesic_Scale(s1,s2,t)
				newTriangles.append(getTriangleFromCan(inter_c, inter_r, inter_t, inter_s))
			mesh.face = np.array(newTriangles)
			mesh.save("/out/path_at_{t}.obj".format())
			print("mesh \'path_at_{}.obj\' saved.".format(t))


			newTriangles = []
			print("Calculating mesh at t=1.25")
			t=1.25
			mesh = geomproc.mesh()
			mesh.vertex = verts1
			newTriangles = []
			for i in range(N):
				c1 = m1.canonical[i]
				c2 = m2.canonical[i]
				r1 = m1.rot[i]
				r2 = m2.rot[i]
				t1 = m1.translation[i]
				t2 = m2.translation[i]
				s1 = m1.scale[i]
				s2 = m2.scale[i] 
				inter_c = geodesic_GA(c1,c2,t)
				inter_r = geodesic_Rot(r1,r2,t)
				inter_t = geodesic_Trans(t1,t2,t)
				inter_s = geodesic_Scale(s1,s2,t)
				newTriangles.append(getTriangleFromCan(inter_c, inter_r, inter_t, inter_s))

			mesh.face = np.array(newTriangles)
			mesh.save("/out/path_at_{t}.obj".format())
			print("mesh \'path_at_{}.obj\' saved.".format(t))

