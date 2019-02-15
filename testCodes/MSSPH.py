import numpy as np
from numpy import newaxis

dx = 0.05
smoothingLength  = 3.0 * dx
smoothingLength2 = 2.0 * dx 

pos = []

def W(i,j,h):
	h = h / 2
	xi     = pos[i]
	xj     = pos[j]
	relpos = xi-xj
	dist   = np.linalg.norm(relpos,2)
	reldir = relpos / (dist + 0.00000000001)
	q = dist / h

	if q > 2.0 :
		return 0
	else:
		foo = (1.0 - 0.5 * q)
		bar = foo * foo
		res = (0.55704230082163367519109317180380 / (h * h)) * bar * bar * (2.0 * q + 1.0)
		return res

def tensorProduct(a,b):
	tens = np.matrix([[0.0,0.0],[0.0,0.0]])
	for _i in range(0,2):
		for _j in range(0,2):
			tens[_i,_j] = a[_i] * b[_j]
	return tens

def gW(i,j,h):
	h = h / 2
	xi     = pos[i]
	xj     = pos[j]
	relpos = xi-xj
	dist   = np.linalg.norm(relpos,2)
	reldir = relpos / (dist + 0.00000000001)
	q = dist / h

	if q > 2.0 :
		return np.array([0.0,0.0])
	else:
		foo = 1.0 - 0.5 * q;
		res = (((7.0/(4.0*np.pi))/(h*h*h))*(-5.0 * q)*foo*foo*foo) * reldir;
		return res

m = 21

for i in range(0,m):
	for j in range(0,m):
		pos.append(np.array([i * dx, j * dx]))

n = len(pos)

f = [0 for _ in range(0,n)]
for i in range(0,n):
	# if (pos[i][0,0] >= 0.3999 and pos[i][0,0] <= 0.60001 and pos[i][1,0] >= 0.3999 and pos[i][1,0] <= 0.60001) :
		# f[i] = 0
	# else:
	f[i] = pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1]

# print("pos[1][0,0] = {}".format(pos[1][0]))
# print(f[i])

kernelSum  = [0 for _ in range(0,n)]
kernelSum2 = [0 for _ in range(0,n)]
vol  = [0 for _ in range(0,n)]
vol2 = [0 for _ in range(0,n)]
L    = [np.matrix([[0.0,0.0],[0.0,0.0]]) for _ in range(0,n)]
L2   = [np.matrix([[0.0,0.0],[0.0,0.0]]) for _ in range(0,n)]
gradf = [np.array([0.0,0.0]) for _ in range(0,n)]
gradf2 = [np.array([0.0,0.0]) for _ in range(0,n)]
gradf_nomralized = [np.array([0.0,0.0]) for _ in range(0,n)]
gradf_nomralized2 = [np.array([0.0,0.0]) for _ in range(0,n)]

laplf_nneighbors = [0 for _ in range(0,n)] 
laplf_nneighbors2 = [0 for _ in range(0,n)] 
laplf_cleary = [0 for _ in range(0,n)]

laplf_nneighbors_normalized = [0 for _ in range(0,n)] 
laplf_nneighbors_normalized_normalized = [0 for _ in range(0,n)] 

laplf_cleary_normalized = [0 for _ in range(0,n)]

for i in range(0,n):
	for j in range(0,n):

		kernelSum[i]  += W(i,j,smoothingLength)
		kernelSum2[i] += W(i,j,smoothingLength2)
		
	vol[i]  = 1./kernelSum[i]
	vol2[i] = 1./kernelSum2[i]
# print(vol)
# print(gradf)
print(np.shape(gradf[1]))
print(np.shape(gW(1,1,0.5)))
for i in range(0,n):
	for j in range(0,n):

		if( np.linalg.norm(pos[i]-pos[j],2) > smoothingLength):
			continue


		gradf [i] += (f[j] - f[i]) * vol [j] * gW(i,j,smoothingLength) 
		gradf2[i] += (f[j] - f[i]) * vol2[j] * gW(i,j,smoothingLength2) 

		L[i]  += tensorProduct((pos[i] - pos[j]), gW(i,j,smoothingLength))  * vol[j]
		L2[i] += tensorProduct((pos[i] - pos[j]), gW(i,j,smoothingLength2)) * vol2[j]

	L [i] = -np.linalg.inv(L [i])
	L2[i] = -np.linalg.inv(L2[i])

	gradf_nomralized[i]  = L[i]  * gradf[i][:,newaxis]
	gradf_nomralized2[i] = L2[i] * gradf2[i][:,newaxis]


print(np.shape(gradf[1]))
print(np.shape(gW(1,1,0.5)))
print(gradf[i])

for i in range(0,n):
	for j in range(0,n):

		relpos = pos[j]-pos[i]
		dist   = np.linalg.norm(relpos,2)
		reldir = relpos / (dist + 0.00000000001)
		# print("asdf")

		# print(gradf2[j]-gradf2[i])
		# print(L[i])
		# print(gW(i,j,smoothingLength))
		# print("fuka")
		# print(np.shape(L[i]),np.shape(gW(i,j,smoothingLength)),np.shape((L[i] * gW(i,j,smoothingLength)[:,newaxis])))
		# print("foobar")
		# print("WTF")
		# print(np.shape( (gradf2[j]-gradf2[i])           ))
		# print(np.shape(  ))


		laplf_cleary[i] 	+= np.dot( 2.0 * reldir * (f[j] - f[i]) / (np.linalg.norm(pos[i]-pos[j],2) + 0.00000000001), gW(i,j,smoothingLength)) * vol[j]
		
		foo = (L[i] * gW(i,j,smoothingLength)[:,newaxis])
		bar = np.array([foo[0,0],foo[1,0]])
		# print("fuck")
		# print(bar)
		# print()
		laplf_nneighbors[i] += np.dot( (gradf2[j]-gradf2[i]), bar ) * vol[j]		
		# print("asdf")
		# print((gradf_nomralized[j].T[0]-gradf_nomralized[i].T[0]))
		# print(gW(i,j,smoothingLength).T)
		# print((np.dot( (gradf_nomralized[j].T[0]-gradf_nomralized[i].T[0]), gW(i,j,smoothingLength).T[0]) * vol[j])[0,0])
		# print(L[i])
		# print(gW(i,j,smoothingLength))
		# print("asdfasdf")
		# print(gradf_nomralized[j].T[0]-gradf_nomralized[i].T[0])
		# print((L[i] * gW(i,j,smoothingLength)).T[0])
		# print("---------")
		# laplf_nneighbors_normalized_normalized[i] += (np.dot( (gradf_nomralized[j].T[0]-gradf_nomralized[i].T[0]).T[0], (L[i] * gW(i,j,smoothingLength)).T[0]) * vol[j])[0,0]
		# laplf_cleary_normalized[i] 	+= np.dot( 2.0 * reldir.T[0] * (f[j] - f[i]) / (np.linalg.norm(pos[i].T[0]-pos[j].T[0],2) + 0.00000000001), gW(i,j,smoothingLength).T[0]) * vol[j]
		

file = open('./test.csv','w')
file.write("x,y,f,gfx,gfx2,gfy,gfy2,laplf_cleary,laplf_nneighbors\n")

for i in range(0,n):
	# print("{},{},{},{},{}".format(pos[i][0],pos[i][1],f[i],gradf[i][0],gradf[i][1]))
	file.write("{},{},{},{},{},{},{},{},{}\n".format(
											   pos[i][0],
											   pos[i][1],
											   f[i],
											   gradf[i][0],
											   gradf2[i][0],
											   gradf[i][1],
											   gradf2[i][1],
											   laplf_cleary[i],
											   laplf_nneighbors[i]))

