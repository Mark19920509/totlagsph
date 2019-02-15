import numpy as np

EPS_SMALL = 1.0E-10

class Particle:
    def __init__(self, idx, x, f):
        self.idx = idx
        self.pos = x
        self.f = f
        self.grad_f = None
        self.vol = None
        self.B = None
        self.L = None
        self.A = None
        self.M = None
        self.N0 = None
        self.N1 = None
        self.N2 = None
        self.neigh = []

    def __str__(self):
        s = "Particle : {}\n".format(idx)
        s+= "vol : " + str(self.vol) + "\n"
        s+= "pos : " + self.pos.__str__() + "\n"
        s+= "f : " + self.f.__str__() + "\n"
        s+= "grad_f : " + self.grad_f.__str__()
        return s

def tensorProduct(a,b) :
    return np.matrix([[a[0] * b[0], a[0] * b[1]],
                      [a[1] * b[0], a[1] * b[1]]])

def W(dist, smoothingLength) :
    H = smoothingLength * 0.5
    q = (dist / H)
    if (q > 2.0) :
 	    return 0;
    else : 
        foo = (1.0 - 0.5 * q);
        bar = foo * foo;
        res = (0.55704230082163367519109317180380 / (h * h)) * bar * bar * (2.0 * q + 1.0)
        return res;

def gW(dist, smoothingLength, direction):
    H = smoothingLength * 0.5;
    q = (dist / H);
    if (q > 2.0) :
        return np.array([0,0])
    else :
        foo = 1.0 - 0.5 * q;
        c = ((0.557042300821633675191)/(H*H*H))*(-5.0 * q)*foo*foo*foo;
        return c * direction;

    


l = 1.0
n = 10

dx = l/float(n)
h = 1.5 * dx

# IDX2D = [(0,0),(1,1),(0,1),(1,0)]
IDX2D = [(0, 0), (0, 1), (1, 1)]
# NEG_DELTA_MN = np.array([-1.0,-1.0,0.0,0.0])
NEG_DELTA_MN = np.array([-1.0, 0.0, -1.0])

particles = []

idx = 0
for i in range(0,n+1):
    for j in range(0,n+1):
        pos = np.array([float(i) * dx, float(j) * dx])
        f = pos[0]
        particles.append(Particle(idx,pos,f))
        idx += 1

# Find Volume
for i in range(0, len(particles)):
    pi = particles[i]
    kernelSum = 0.0
    for j in range(0, len(particles)):
        pj = particles[j]

        dist_ij = np.linalg.norm((pi.pos - pj.pos), 2)
        if dist_ij > h :
            continue        

        # e_ij = (pi.pos - pj.pos) / dist_ij
        W_ij = W(dist_ij, h)
        # gW_ij = gW(dist,h,e_ij)

        kernelSum += W_ij        
    pi.vol = 1.0 / kernelSum

# Find Renorm Matrix (1st derivative)
for i in range(0, len(particles)):
    pi = particles[i]
    vol_i = pi.vol
    f_i = pi.f
    B_i = np.matrix([[0,0],[0,0]])
    
    grad_f_i = np.array([0,0])
    for j in range(0, len(particles)):
        pj = particles[j]

        dist_ij = np.linalg.norm((pi.pos - pj.pos), 2)
        if ( dist_ij > h or pi == pj):
            continue        
        vol_j = pj.vol
        f_j = pj.f

        r_ij = (pi.pos - pj.pos)
        e_ij = (pi.pos - pj.pos) / dist_ij
        W_ij = W(dist_ij, h)
        gW_ij = gW(dist_ij,h,e_ij)
        
        B_i = B_i - vol_j * tensorProduct(r_ij, gW_ij)
        grad_f_i = grad_f_i + (vol_j * (f_j - f_i) * gW_ij)

    pi.B = B_i = np.linalg.inv(B_i)
    pi.grad_f = B_i.dot(grad_f_i.transpose())





# Find Renorm Matrix (2nd derivative)
for i in range(0, len(particles)):
    pi = particles[i]
    vol_i = pi.vol
    f_i = pi.f
    grad_f_i = pi.grad_f
    B_i = pi.B

    A_i = np.einsum('m,n,q->mnq',np.array([0,0]),np.array([0,0]),np.array([0,0]))
    for j in range(0, len(particles)):
        pj = particles[j]

        dist_ij = np.linalg.norm((pi.pos - pj.pos), 2)
        if ( dist_ij > h or pi == pj):
            continue        
        vol_j = pj.vol
        f_j = pj.f

        r_ij = (pi.pos - pj.pos)
        e_ij = (pi.pos - pj.pos) / dist_ij
        W_ij = W(dist_ij, h)
        gW_ij = gW(dist_ij,h,e_ij)

        G_i_ok_gW_k = np.einsum('ok,k->o',B_i,gW_ij)
        A_i = A_i + vol_j * np.einsum('m,n,o->mno',r_ij,r_ij,G_i_ok_gW_k)
        # M_i_mnq = M_i_mnq + vol_j * np.einsum('m,n,q->mnq',r_ij,r_ij,gW_ij)

    # pi.A = np.einsum('kq,mnq->kmn',B_i,M_i_mnq)    
    pi.A = A_i



# Find Renorm Matrix (2nd derivative)
# for i in range(0, len(particles)):
#     pi = particles[i]
#     vol_i = pi.vol
#     f_i = pi.f
#     grad_f_i = pi.grad_f
#     B_i = pi.B

#     N0_i = np.einsum('m,n->mn',np.array([0,0]),np.array([0,0]))
#     N1_i = np.einsum('m,n->mn',np.array([0,0]),np.array([0,0]))
#     N2_i = np.einsum('m,n->mn',np.array([0,0]),np.array([0,0]))

#     for j in range(0, len(particles)):
#         pj = particles[j]

#         dist_ij = np.linalg.norm((pi.pos - pj.pos), 2)
#         if ( dist_ij > h or pi == pj):
#             continue        
#         vol_j = pj.vol
#         f_j = pj.f

#         r_ij = (pi.pos - pj.pos)
#         e_ij = (pi.pos - pj.pos) / dist_ij
#         W_ij = W(dist_ij, h)
#         gW_ij = gW(dist_ij,h,e_ij)

#         N0_i = N1_i + vol_j * np.einsum('m,n->mn',r_ij,gW_ij)
#         N1_i = N1_i + vol_j * np.einsum('m,n->mn',e_ij,gW_ij)
#         N2_i = N2_i + vol_j * np.einsum('m,n->mn',r_ij,r_ij * np.norm(gW_ij, 2))

#     pi.A = np.einsum('kq,mnq->kmn',B_i,M_i_mnq)    
#     pi.N0 = N0_i
#     pi.N1 = N1_i
#     pi.N2 = N2_i


# for i in range(0, len(particles)):
#     pi = particles[i]
#     vol_i = pi.vol
#     f_i = pi.f
#     grad_f_i = pi.grad_f
#     B_i = pi.B

#     M_i = np.einsum('m,n,o->mno',np.array([0,0]),np.array([0,0]),np.array([0,0]))
#     for j in range(0, len(particles)):
#         pj = particles[j]

#         dist_ij = np.linalg.norm((pi.pos - pj.pos), 2)
#         if ( dist_ij > h or pi == pj):
#             continue        
#         vol_j = pj.vol
#         f_j = pj.f

#         r_ij = (pi.pos - pj.pos)
#         e_ij = (pi.pos - pj.pos) / dist_ij
#         W_ij = W(dist_ij, h)
#         gW_ij = gW(dist_ij,h,e_ij)
        
#         M_i = M_i + vol_j * np.einsum('m,n->mn',r_ij,gW_ij) + 


# Find Renorm Matrix (continued)

# for i in range(0, len(particles)):
#     pi = particles[i]

#     vol_i = pi.vol
#     f_i = pi.f
#     grad_f_i = pi.grad_f
#     B_i = pi.B
#     A_i = pi.A
    
#     # N_i = np.einsum('m,n,o,p->mnop',np.array([0,0]),np.array([0,0]),np.array([0,0]),np.array([0,0]))
#     nneigh = 0
#     for j in range(0, len(particles)):
#         pj = particles[j]

#         dist_ij = np.linalg.norm((pi.pos - pj.pos), 2)
#         if ( dist_ij > h or pi == pj):
#             continue        
#         nneigh += 1
#         vol_j = pj.vol
#         f_j = pj.f

#         r_ij = (pi.pos - pj.pos)
#         e_ij = (pi.pos - pj.pos) / dist_ij
#         W_ij = W(dist_ij, h)
#         gW_ij = gW(dist_ij,h,e_ij)

#         A_i_kmn_e_ij_k = np.einsum('kmn,k->mn',A_i,e_ij)
#         r_ij_e_ij_mn = np.einsum('m,n->mn',r_ij,e_ij)             
#         e_ij_gW_ij_op = np.einsum('o,p->op',e_ij,gW_ij)   
#         _N_i = _N_i + vol_j * np.einsum('mn,op->mnop',( A_i_kmn_e_ij_k + r_ij_e_ij_mn ), e_ij_gW_ij_op)

#     N_i = np.einsum('I,J->IJ',np.array([0,0,0,0]),np.array([0,0,0,0]))
#     for I in range(len(IDX2D)):
#         m = IDX2D[I][0]
#         n = IDX2D[I][1]

#         for J in range(len(IDX2D)):
#             o = IDX2D[J][0]
#             p = IDX2D[J][1]

#             N_i[I,J] = _N_i[m,n,o,p]
#     print(nneigh)
    # print(np.linalg.inv(N_i).dot(NEG_DELTA_MN))
    




for i in range(0, len(particles)):
    pi = particles[i]

    vol_i = pi.vol
    f_i = pi.f
    grad_f_i = pi.grad_f
    B_i = pi.B
    A_i = pi.A
    
    N_i = np.einsum('i,j->ij',np.array([0,0,0]),np.array([0,0,0]))
    print(pi)
    for j in range(0, len(particles)):
        pj = particles[j]

        dist_ij = np.linalg.norm((pi.pos - pj.pos), 2)
        if ( dist_ij > h or pi == pj):
            continue        
        vol_j = pj.vol
        f_j = pj.f

        r_ij = (pi.pos - pj.pos)
        e_ij = (pi.pos - pj.pos) / dist_ij
        W_ij = W(dist_ij, h)
        gW_ij = gW(dist_ij,h,e_ij)

        A_i_kmn_e_ij_k = np.einsum('kmn,k->mn',A_i,e_ij)
        r_ij_e_ij_mn = np.einsum('m,n->mn',r_ij,e_ij)             
        e_ij_gW_ij_op = np.einsum('o,p->op',e_ij,gW_ij)   
        
        # (0,0)(0,1)(1,1)
        for I in range(0,len(IDX2D)):
            m = IDX2D[I][0]
            n = IDX2D[I][1]

            for J in range(0,len(IDX2D)):
                o = IDX2D[J][0]
                p = IDX2D[J][1]

                if o == p :
                    N_i[I,J] = N_i[I,J] + vol_j * ( A_i_kmn_e_ij_k[m,n] + r_ij_e_ij_mn[m,n] ) * (e_ij_gW_ij_op[o,p])
                else :
                    N_i[I,J] = N_i[I,J] + vol_j * ( A_i_kmn_e_ij_k[m,n] + r_ij_e_ij_mn[m,n] ) * (e_ij_gW_ij_op[o,p] + e_ij_gW_ij_op[p,o])

    print(N_i)
    # print(np.linalg.inv(N_i).dot(NEG_DELTA_MN))
