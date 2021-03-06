import  numpy          as      np
from    .configuration  import  *
from    .coordParser    import  wrapAngles, calcNorm
from    math           import  sqrt
from    sys  import  stderr
from    copy import deepcopy
from    .fastmath import fastCross

def c6c12_to_sigmaepsilon(c6, c12):
    return ((c12/c6)**(1.0/6), 0.25*(c6**2/c12))

def sigmaepsilon_to_c6c12(sigma, epsilon):
    return (4*epsilon*(sigma**6), 4*epsilon*(sigma**12))

class atomTerms (object):

    def __init__ (self, size):
        self.size = size
        self.c6 = np.zeros(size)
        self.c12 = np.zeros(size)
        self.cs6 = np.zeros(size)
        self.cs12 = np.zeros(size)
        self.q = np.zeros(size)

    def setMember (self, i, c6, c12, cs6, cs12, q):
        self.c6[i] = c6
        self.c12[i] = c12
        self.cs6[i] = cs6
        self.cs12[i] = cs12
        self.q[i] = q

    def mixLJ (self, i, j, type, mixtype='geometric'):
        if (mixtype == 'geometric'):
            if (type == 1):
                return (sqrt(self.c6[i]*self.c6[j]), sqrt(self.c12[i]*self.c12[j]))
            if (type == 2):
                return (sqrt(self.cs6[i]*self.cs6[j]), sqrt(self.cs12[i]*self.cs12[j]))
        if (mixtype == 'arithmetic'):
            if (type == 1):
                sigma_i, epsilon_i = c6c12_to_sigmaepsilon(self.c6[i], self.c12[i])
                sigma_j, epsilon_j = c6c12_to_sigmaepsilon(self.c6[j], self.c12[j])
            if (type == 2):
                sigma_i, epsilon_i = c6c12_to_sigmaepsilon(self.cs6[i], self.cs12[i])
                sigma_j, epsilon_j = c6c12_to_sigmaepsilon(self.cs6[j], self.cs12[j])
            sigma = 0.50 * (sigma_i + sigma_j)
            epsilon = np.sqrt(epsilon_i * epsilon_j)
            c6, c12 = sigmaepsilon_to_c6c12(sigma, epsilon)
            return (c6, c12)

    def mixQ (self, i, j):
        return self.q[i] * self.q[j]

class G96bondTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(size, dtype=np.uint8)
        self.aj = np.zeros(size, dtype=np.uint8)
        self.k  = np.zeros(size)
        self.l  = np.zeros(size)

    def setMember (self, i, ai, aj, k, l):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.k[i] = k
        self.l[i] = l

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0
        difference2  = conf.getDistances2(self.ai, self.aj)
        difference2 -= self.l * self.l
        energies = 0.25 * self.k * difference2 * difference2
        rijs = conf.getDisplacements (self.ai, self.aj)
        f_i  = np.transpose(-self.k * difference2 * np.transpose(rijs))
        for ib in range(self.size):
            i = self.ai[ib]
            j = self.aj[ib]
            forceConf[i,:] +=  f_i[ib,:]
            forceConf[j,:] += -f_i[ib,:]
        return energies

class harmonicBondTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(size, dtype=np.uint8)
        self.aj = np.zeros(size, dtype=np.uint8)
        self.k  = np.zeros(size)
        self.l  = np.zeros(size)

    def setMember (self, i, ai, aj, k, l):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.k[i] = k
        self.l[i] = l

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0
        difference  = conf.getDistances(self.ai, self.aj)
        difference2 = conf.getDistances2(self.ai, self.aj)
        rijs_m = deepcopy(difference)
        difference -= self.l
        energies = 0.5 * self.k * (difference ** 2)
        rijs = conf.getDisplacements (self.ai, self.aj)
        f_i  = np.transpose(-self.k * difference * np.transpose(rijs) / rijs_m)
        for ib in range(self.size):
            i = self.ai[ib]
            j = self.aj[ib]
            forceConf[i,:] +=  f_i[ib,:]
            forceConf[j,:] += -f_i[ib,:]
        return energies
    
class G96angleTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(size, dtype=np.uint8)
        self.aj = np.zeros(size, dtype=np.uint8)
        self.ak = np.zeros(size, dtype=np.uint8)
        self.k  = np.zeros(size)
        self.theta = np.zeros(size)

    def setMember (self, i, ai, aj, ak, k, theta):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.ak[i] = ak - 1
        self.k[i]  = k
        self.theta[i] = theta

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0.0
        cos_theta0 = np.cos(np.radians(self.theta))
        cos_theta = conf.getCosines (self.ai, self.aj, self.ak)
        diff_cos = cos_theta - cos_theta0
        energies = 0.5 * self.k * diff_cos * diff_cos
        # forces
        mod_rij  = conf.getDistances(self.ai, self.aj)
        mod_rkj  = conf.getDistances(self.ak, self.aj)
        rij_norm = np.transpose(conf.getDisplacements(self.ai, self.aj)) / mod_rij
        rkj_norm = np.transpose(conf.getDisplacements(self.ak, self.aj)) / mod_rkj
        f_i = np.transpose( -self.k * diff_cos * (rkj_norm - cos_theta*(rij_norm)) / mod_rij )
        f_k = np.transpose( -self.k * diff_cos * (rij_norm - cos_theta*(rkj_norm)) / mod_rkj )
        for ia in range(self.size):
            forceConf[self.ai[ia],:] += f_i[ia,:]
            forceConf[self.ak[ia],:] += f_k[ia,:]
            forceConf[self.aj[ia],:] += -(f_i[ia,:] + f_k[ia,:])
        return energies

class harmonicAngleTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(size, dtype=np.uint8)
        self.aj = np.zeros(size, dtype=np.uint8)
        self.ak = np.zeros(size, dtype=np.uint8)
        self.k  = np.zeros(size)
        self.theta = np.zeros(size)

    def setMember (self, i, ai, aj, ak, k, theta):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.ak[i] = ak - 1
        self.k[i]  = k
        self.theta[i] = theta

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0.0
        theta0 = np.radians(self.theta)
        theta  = np.radians(conf.getAngles(self.ai, self.aj, self.ak))
        cos_theta    = conf.getCosines(self.ai, self.aj, self.ak)
        sin_theta    = conf.getSines(self.ai, self.aj, self.ak)
        diff_theta = (theta - theta0)
        energies = 0.5 * self.k * diff_theta * diff_theta
        # forces
        mod_rij  = conf.getDistances(self.ai, self.aj)
        mod_rkj  = conf.getDistances(self.ak, self.aj)
        rij_norm = np.transpose(conf.getDisplacements(self.ai, self.aj)) / mod_rij
        rkj_norm = np.transpose(conf.getDisplacements(self.ak, self.aj)) / mod_rkj
        f_i = np.transpose( self.k * (diff_theta/sin_theta) * (rkj_norm - cos_theta*(rij_norm)) / mod_rij )
        f_k = np.transpose( self.k * (diff_theta/sin_theta) * (rij_norm - cos_theta*(rkj_norm)) / mod_rkj )
        for ia in range(self.size):
            forceConf[self.ai[ia],:] += f_i[ia,:]
            forceConf[self.ak[ia],:] += f_k[ia,:]
            forceConf[self.aj[ia],:] += -(f_i[ia,:] + f_k[ia,:])
        return energies

class RBAngleTerms (object):
    """Composition of harmonic angle and harmonic bond."""
    def __init__ (self, size):
        self.harmonicAngleTerms = harmonicAngleTerms(size)
        self.harmonicBondTerms  = harmonicBondTerms(size)

    def setMember (self, i, ai, aj, ak, k, theta, kub, r13):
        self.harmonicAngleTerms.setMember(i, ai, aj, ak, k, theta)
        self.harmonicBondTerms.setMember(i, ai, ak, kub, r13)

    def calcForConf (self, conf, forceConf):
        ea = self.harmonicAngleTerms.calcForConf(conf, forceConf)
        eb = self.harmonicBondTerms.calcForConf(conf, forceConf)
        return (ea + eb)

class dihedralTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(self.size, dtype=np.uint8)
        self.aj = np.zeros(self.size, dtype=np.uint8)
        self.ak = np.zeros(self.size, dtype=np.uint8)
        self.al = np.zeros(self.size, dtype=np.uint8)
        self.phi = np.zeros(self.size)
        self.k = np.zeros(self.size)
        self.m = np.zeros(self.size)

    def getType(self):
        return "standard"

    def verifyPhase (self, i):
        if (self.phi[i] != 0) and (self.phi[i] != 180.00):
            raise RuntimeError ("For symmetric dihedral, phase must be 0 or 180.0 deg; it is {:.2f}".format(self.phi[i]))

    def setMember (self, i, ai, aj, ak, al, phi, k, m):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.ak[i] = ak - 1
        self.al[i] = al - 1
        self.phi[i] = phi
        self.k[i] = k
        self.m[i] = m
        self.verifyPhase(i)

    def duplicateDihedral (self, i):
        self.ai = np.append(self.ai, np.array([self.ai[i]]))
        self.aj = np.append(self.aj, np.array([self.aj[i]]))
        self.ak = np.append(self.ak, np.array([self.ak[i]]))
        self.al = np.append(self.al, np.array([self.al[i]]))
        self.phi = np.append(self.phi, np.array([self.phi[i]]))
        self.k = np.append(self.k, np.array([self.k[i]]))
        self.m = np.append(self.m, np.array([self.m[i]]))
        self.size += 1
        return self.size - 1

    def copy (self):
        newTerm = dihedralTerm (self.size)
        for i in range(self.size):
            newTerm.setMember(i, self.ai[i], self.aj[i], self.ak[i], self.al[i],
                    self.phi[i], self.k[i], self.m[i])
        return newTerm

    def setParameters (self, i, phi=None, k=None, m=None):
        if phi is not None:
            self.phi[i] = phi
            self.verifyPhase(i)
        if k is not None:
            self.k[i] = k
        if m is not None:
            self.m[i] = m

    # this auxiliary function returns the value of 
    # d(cos(m*phi))/d(cos(phi))
    @staticmethod
    def cosineDerivatives (cos_phi, m):
        c = cos_phi
        if (m == 0):
            return 0
        elif (m == 1):
            return 1
        elif (m == 2):
            return 4 * c
        elif (m == 3):
            return 12 * (c*c) - 3
        elif (m == 4):
            return 32 * (c*c*c) - 16 * c
        elif (m == 5):
            return 80 * (c*c*c*c) - 60 * (c*c) + 5
        elif (m == 6):
            return 192 * (c*c*c*c*c) - 192 * (c*c*c) + 36 * c
        else:
            raise RuntimeError ("Only multiplicities 0-6 are supported.")

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0.0
        phi = np.radians(conf.getDihedrals(self.ai, self.aj, self.ak, self.al))
        cos_phi_0 = np.cos(np.radians(self.phi))
        energies = self.k * (1 + cos_phi_0 * np.cos(self.m * phi))
        # force
        rij = conf.getDisplacements(self.ai, self.aj)
        rkj = conf.getDisplacements(self.ak, self.aj)
        rkl = conf.getDisplacements(self.ak, self.al)
        mod_rkj = conf.getDistances(self.ak, self.aj)
        mod_rkj2 = mod_rkj * mod_rkj
        dot_ij_kj = np.einsum('ij,ij->i', rij, rkj)
        dot_kl_kj = np.einsum('ij,ij->i', rkl, rkj)
        # auxiliary vectors
        rim =  rij - np.transpose(dot_ij_kj * np.transpose(rkj) / (mod_rkj2))
        rln = -rkl + np.transpose(dot_kl_kj * np.transpose(rkj) / (mod_rkj2))
        mod_rim = np.sqrt( np.einsum('ij,ij->i', rim, rim) )
        mod_rln = np.sqrt( np.einsum('ij,ij->i', rln, rln) )
        rim_norm = np.transpose(rim) / mod_rim # ALREADY TRANSPOSED
        rln_norm = np.transpose(rln) / mod_rln # ALREADY TRANSPOSED
        #
        cos_phi = np.cos(phi)
        # cosine derivatives
        dcos = np.array([dihedralTerms.cosineDerivatives(cos_phi[i], self.m[i]) for
            i in range(self.size)])
        # forces
        f_i = -self.k*cos_phi_0*dcos*(rln_norm - cos_phi*rim_norm) / (mod_rim)
        f_l = -self.k*cos_phi_0*dcos*(rim_norm - cos_phi*rln_norm) / (mod_rln)
        f_j = ( (dot_ij_kj/(mod_rkj2)) - 1 ) * f_i - (dot_kl_kj/(mod_rkj2)) * f_l
        f_k = -(f_i + f_l + f_j)
        for idih in range(self.size):
            forceConf[self.ai[idih],:] += np.transpose(f_i)[idih,:]
            forceConf[self.aj[idih],:] += np.transpose(f_j)[idih,:]
            forceConf[self.ak[idih],:] += np.transpose(f_k)[idih,:]
            forceConf[self.al[idih],:] += np.transpose(f_l)[idih,:]
        return energies

class optDihedralTerms(object):
    """
    Special optDihedralTerm for optimized dihedrals.

    The difference between this one and the usual dihedralTerm is
    that, for performance reasons, this term includes the sum of all
    multiplicities, while the usual term refers to only one
    multiplicity. 
    """

    def __init__(self, size):
        self.size   = size
        self.ai     = np.zeros(self.size, dtype=np.uint8)
        self.aj     = np.zeros(self.size, dtype=np.uint8)
        self.ak     = np.zeros(self.size, dtype=np.uint8)
        self.al     = np.zeros(self.size, dtype=np.uint8)
        self.k      = np.zeros((6, self.size))
        self.phi    = np.zeros((6, self.size))
        self.m      = np.array([1,2,3,4,5,6], dtype=np.uint8)

    def verifyPhase (self):
        for phi in self.phi.flatten():
            if (phi != 0) and (phi != 180.00):
                raise RuntimeError ("For symmetric dihedral, phase must be 0 or 180.0 deg; it is {:.2f}".format(self.phi[i]))

    def getType(self):
        return "standard"

    def setMember(self, i, ai, aj, ak, al):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.ak[i] = ak - 1
        self.al[i] = al - 1

    def setParameters (self, m, phi, k):
        if m < 1:
            raise Exception
        if k is not None:
            self.k[m-1,:] = k
        if phi is not None:
            self.phi[m-1,:] = phi
        self.verifyPhase()

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0.0
        phi       = np.radians(conf.getDihedrals(self.ai, self.aj, self.ak, self.al))
        cos_phi_0 = np.cos(np.radians(self.phi))
        energies  = np.einsum('ij,ik->k', self.k, (1 + cos_phi_0 * np.cos(np.einsum('i,j->ij', self.m, phi))))
        # force
        rij = conf.getDisplacements(self.ai, self.aj)
        rkj = conf.getDisplacements(self.ak, self.aj)
        rkl = conf.getDisplacements(self.ak, self.al)
        mod_rkj = conf.getDistances(self.ak, self.aj)
        mod_rkj2 = mod_rkj * mod_rkj
        dot_ij_kj = np.einsum('ij,ij->i', rij, rkj)
        dot_kl_kj = np.einsum('ij,ij->i', rkl, rkj)
        # auxiliary vectors
        rim =  rij - np.transpose(dot_ij_kj * np.transpose(rkj) / (mod_rkj2))
        rln = -rkl + np.transpose(dot_kl_kj * np.transpose(rkj) / (mod_rkj2))
        mod_rim = np.sqrt( np.einsum('ij,ij->i', rim, rim) )
        mod_rln = np.sqrt( np.einsum('ij,ij->i', rln, rln) )
        rim_norm = np.transpose(rim) / mod_rim # ALREADY TRANSPOSED
        rln_norm = np.transpose(rln) / mod_rln # ALREADY TRANSPOSED
        #
        cos_phi = np.cos(phi)
        # cosine derivatives
        dcos = np.array([[dihedralTerms.cosineDerivatives(cos_phi[j], m) for j in range(self.size)] for m in self.m])
        # forces


        f_i = np.einsum('ij,kj->kj', -self.k*cos_phi_0*dcos, (rln_norm - cos_phi*rim_norm) / (mod_rim))
        f_l = np.einsum('ij,kj->kj', -self.k*cos_phi_0*dcos, (rim_norm - cos_phi*rln_norm) / (mod_rln))
        f_j = ( (dot_ij_kj/(mod_rkj2)) - 1 ) * f_i - (dot_kl_kj/(mod_rkj2)) * f_l
        f_k = -(f_i + f_l + f_j)
        for idih in range(len(phi)):
            forceConf[self.ai[idih],:] += np.transpose(f_i)[idih,:]
            forceConf[self.aj[idih],:] += np.transpose(f_j)[idih,:]
            forceConf[self.ak[idih],:] += np.transpose(f_k)[idih,:]
            forceConf[self.al[idih],:] += np.transpose(f_l)[idih,:]
        return energies

class RyckaertBellemansDihedralTerms(object):

    def __init__(self, size):
        self.size = size
        self.ai = np.zeros(self.size, dtype=np.uint8)
        self.aj = np.zeros(self.size, dtype=np.uint8)
        self.ak = np.zeros(self.size, dtype=np.uint8)
        self.al = np.zeros(self.size, dtype=np.uint8)
        self.c0 = np.zeros(self.size)
        self.c1 = np.zeros(self.size)
        self.c2 = np.zeros(self.size)
        self.c3 = np.zeros(self.size)
        self.c4 = np.zeros(self.size)
        self.c5 = np.zeros(self.size)

    def getType(self):
        return "ryckaert"

    def setMember(self, i, ai, aj, ak, al, c0, c1, c2, c3, c4, c5):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.ak[i] = ak - 1
        self.al[i] = al - 1
        self.c0[i] = c0
        self.c1[i] = c1
        self.c2[i] = c2
        self.c3[i] = c3
        self.c4[i] = c4
        self.c5[i] = c5

    def setParameters (self, i, j, c):
        if (j < 0):
            raise Exception()
        if c is not None:
            if j == 0:
                self.c0[i] = c
            if j == 1:
                self.c1[i] = c
            if j == 2:
                self.c2[i] = c
            if j == 3:
                self.c3[i] = c
            if j == 4:
                self.c4[i] = c
            if j == 5:
                self.c5[i] = c

    def calcForConf(self, conf, forceConf):
        if (self.size == 0):
            return 0
        
        phi = np.radians(conf.getDihedrals(self.ai, self.aj, self.ak, self.al))
        cos = np.cos(phi)
        sin = np.sin(phi)

        energies = self.c0 - self.c1 * cos + self.c2 * (cos**2) - self.c3 * (cos**3) + self.c4 * (cos**4) - self.c5 * (cos**5)
        # force
        force_pre_factor = sin * (-self.c1 + 2*self.c2*cos -3*self.c3*(cos**2) + 4*self.c4*(cos**3) -5*self.c5*(cos**4))
        rij = conf.getDisplacements(self.ai, self.aj)
        rkj = conf.getDisplacements(self.ak, self.aj)
        rkl = conf.getDisplacements(self.ak, self.al)
        # cross products - works for arrays of vectors too!
        rmj = fastCross(rij, rkj)
        rnk = fastCross(rkj, rkl)
        # norms
        mod_rmj = np.linalg.norm(rmj, axis=1)
        mod_rnk = np.linalg.norm(rnk, axis=1)
        mod_rkj = np.linalg.norm(rkj, axis=1)
        mod_rkj2 = mod_rkj * mod_rkj
        # forces
        f_i = +force_pre_factor * (mod_rkj / ((mod_rmj * mod_rmj))) * np.transpose(rmj) # transpose!
        f_l = -force_pre_factor * (mod_rkj / ((mod_rnk * mod_rnk))) * np.transpose(rnk) # transpose!
        f_j = ( (np.einsum('ij,ij->i',rij,rkj)/(mod_rkj2)) - 1 ) * f_i - (np.einsum('ij,ij->i',rkl,rkj)/(mod_rkj2)) * f_l # NO transpose
        f_k = -(f_i + f_l + f_j) # no transpose
        for ir in range(self.size):
            forceConf[self.ai[ir],:] += np.transpose(f_i)[ir,:]
            forceConf[self.aj[ir],:] += np.transpose(f_j)[ir,:]
            forceConf[self.ak[ir],:] += np.transpose(f_k)[ir,:]
            forceConf[self.al[ir],:] += np.transpose(f_l)[ir,:]
        return energies

class FourierDihedralTerms(RyckaertBellemansDihedralTerms):
    def setMember(self, i, ai, aj, ak, al, f1, f2, f3, f4):
        c0 = f2 + 0.50 * (f1 + f3)
        c1 = 0.50 * (-f1 + 3*f3)
        c2 = -f2 + 4*f4
        c3 = -2*f3
        c4 = -4*f4
        c5 = 0
        super().setMember(i, ai, aj, ak, al, c0, c1, c2, c3, c4, c5)

class improperTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(self.size, dtype=np.uint8)
        self.aj = np.zeros(self.size, dtype=np.uint8)
        self.ak = np.zeros(self.size, dtype=np.uint8)
        self.al = np.zeros(self.size, dtype=np.uint8)
        self.phi = np.zeros(self.size)
        self.k = np.zeros(self.size)

    def setMember (self, i, ai, aj, ak, al, k, phi):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.ak[i] = ak - 1
        self.al[i] = al - 1
        self.phi[i] = phi
        self.k[i] = k

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0
        angDiff = wrapAngles(conf.getImpropers(self.ai, self.aj, self.ak, self.al) - self.phi)
        energies = 0.50 * self.k * angDiff * angDiff
        # force
        rij = conf.getDisplacements(self.ai, self.aj)
        rkj = conf.getDisplacements(self.ak, self.aj)
        rkl = conf.getDisplacements(self.ak, self.al)
        # cross products - works for arrays of vectors too!
        rmj = fastCross(rij, rkj)
        rnk = fastCross(rkj, rkl)
        # norms
        mod_rmj = np.linalg.norm(rmj, axis=1)
        mod_rnk = np.linalg.norm(rnk, axis=1)
        mod_rkj = np.linalg.norm(rkj, axis=1)
        mod_rkj2 = mod_rkj * mod_rkj
        # forces
        f_i = -self.k * angDiff * (mod_rkj / ((mod_rmj * mod_rmj))) * np.transpose(rmj) # transpose!
        f_l = +self.k * angDiff * (mod_rkj / ((mod_rnk * mod_rnk))) * np.transpose(rnk) # transpose!
        f_j = ( (np.einsum('ij,ij->i',rij,rkj)/(mod_rkj2)) - 1 ) * f_i - (np.einsum('ij,ij->i',rkl,rkj)/(mod_rkj2)) * f_l # NO transpose
        f_k = -(f_i + f_l + f_j) # no transpose
        for ir in range(self.size):
            forceConf[self.ai[ir],:] += np.transpose(f_i)[ir,:]
            forceConf[self.aj[ir],:] += np.transpose(f_j)[ir,:]
            forceConf[self.ak[ir],:] += np.transpose(f_k)[ir,:]
            forceConf[self.al[ir],:] += np.transpose(f_l)[ir,:]
        return energies

class periodicImproperTerms (dihedralTerms): pass

class dihedralRestraintTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(self.size, dtype=np.uint8)
        self.aj = np.zeros(self.size, dtype=np.uint8)
        self.ak = np.zeros(self.size, dtype=np.uint8)
        self.al = np.zeros(self.size, dtype=np.uint8)
        self.phi = np.zeros(self.size)
        self.k = np.zeros(self.size)

    def setMember (self, i, ai, aj, ak, al, phi, k):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.ak[i] = ak - 1
        self.al[i] = al - 1
        self.phi[i] = phi
        self.k[i] = k

    def popMember (self):
        self.size -= 1
        self.ai = np.delete(self.ai, -1)
        self.aj = np.delete(self.aj, -1)
        self.ak = np.delete(self.ak, -1)
        self.al = np.delete(self.al, -1)
        self.phi = np.delete(self.phi, -1)
        self.k = np.delete(self.k, -1)

    def pushMember (self, ai, aj, ak, al, phi, k):
        self.size += 1
        self.ai = np.append(self.ai, 0)
        self.aj = np.append(self.aj, 0)
        self.ak = np.append(self.ak, 0)
        self.al = np.append(self.al, 0)
        self.phi = np.append(self.phi, 0)
        self.k = np.append(self.k, 0)
        self.setMember(self.size - 1, ai, aj, ak, al, phi, k)

    def print(self):
        print("{} restraints:".format(self.size))
        for i in range(self.size):
            print("ai = {}, aj = {}, ak = {}, al = {}, phi = {}, k = {}".format(
                self.ai[i], self.aj[i], self.ak[i], self.al[i], self.phi[i], self.k[i]))

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0
        angDiff = wrapAngles(conf.getImpropers(self.ai, self.aj, self.ak, self.al) - self.phi)
        energies = 0.50 * self.k * angDiff * angDiff
        # force
        rij = conf.getDisplacements(self.ai, self.aj)
        rkj = conf.getDisplacements(self.ak, self.aj)
        rkl = conf.getDisplacements(self.ak, self.al)
        # cross products - works for arrays of vectors too!
        rmj = fastCross(rij, rkj)
        rnk = fastCross(rkj, rkl)
        # norms
        mod_rmj = np.linalg.norm(rmj, axis=1)
        mod_rnk = np.linalg.norm(rnk, axis=1)
        mod_rkj = np.linalg.norm(rkj, axis=1)
        mod_rkj2 = mod_rkj * mod_rkj
        # forces
        f_i = -self.k * angDiff * (mod_rkj / ((mod_rmj * mod_rmj))) * np.transpose(rmj) # transpose!
        f_l = +self.k * angDiff * (mod_rkj / ((mod_rnk * mod_rnk))) * np.transpose(rnk) # transpose!
        f_j = ( (np.einsum('ij,ij->i',rij,rkj)/(mod_rkj2)) - 1 ) * f_i - (np.einsum('ij,ij->i',rkl,rkj)/(mod_rkj2)) * f_l # NO transpose
        f_k = -(f_i + f_l + f_j) # no transpose
        for ir in range(self.size):
            forceConf[self.ai[ir],:] += np.transpose(f_i)[ir,:]
            forceConf[self.aj[ir],:] += np.transpose(f_j)[ir,:]
            forceConf[self.ak[ir],:] += np.transpose(f_k)[ir,:]
            forceConf[self.al[ir],:] += np.transpose(f_l)[ir,:]
        return energies

class LJTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(self.size, dtype=np.uint8)
        self.aj = np.zeros(self.size, dtype=np.uint8)
        self.type = np.zeros(self.size, dtype=np.uint8)
        self.c6 = np.zeros(self.size)
        self.c12 = np.zeros(self.size)

    def setMember (self, i, ai, aj, type, c6, c12):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.type[i] = type
        self.c6[i] = c6
        self.c12[i] = c12

    def print(self):
        for i in range(self.size):
            print("LJ interaction, type = %d, %d-%d, c6 = %e, c12 = %e" % 
                    (self.type[i], self.ai[i]+1, self.aj[i]+1, self.c6[i], self.c12[i]))

    def setParameters(self, i, c6=None, c12=None):
        if c6 is not None:
            self.c6[i] = c6
        if c12 is not None:
            self.c12[i] = c12

    def setFromAtoms (self, atomterms, mixtype='geometric'):
        for i in range(self.size):
            pars = atomterms.mixLJ(self.ai[i], self.aj[i], self.type[i], mixtype)
            self.c6[i] = pars[0]
            self.c12[i] = pars[1]

    def calcForConf (self, conf, forceConf, debug=False):
        if (self.size == 0):
            return 0
        dist = conf.getDistances (self.ai, self.aj)
        invDist6 = dist ** (-6)
        energies = -self.c6 * ( invDist6 ) + self.c12 * ( invDist6 * invDist6 )

        if (debug):
            sum_14 = 0
            sum_SR = 0
            sum_totl = 0
            print("------------------------------------------------------------------------\n", file=stderr)
            print(" LJ TERMS CALCULATION\n", file=stderr)
            print("------------------------------------------------------------------------\n", file=stderr)
            print("NTERMS = %d\n" % self.size, file=stderr)
            for i in range(self.size):
                disp = -self.c6[i]*invDist6[i]
                repl = +self.c12[i]*(invDist6[i]*invDist6[i])
                totl = disp + repl
                sum_totl += totl
                if (self.type[i] == 1):
                    sum_SR += totl
                elif (self.type[i] == 2):
                    sum_14 += totl
                else:
                    raise ValueError()
                print("i=%d, j=%d, type=%d" % (self.ai[i] + 1, self.aj[i] + 1, self.type[i]), file=stderr)
                print("    rij[%d] =%10.4f" % (i, dist[i]), file=stderr)
                print("  rij-6[%d] =%10.4f" % (i, invDist6[i]), file=stderr)
                print("   c6ij[%d] =%18.7e" % (i, self.c6[i]), file=stderr)
                print("  c12ij[%d] =%18.7e" % (i, self.c12[i]), file=stderr)
                print("   disp[%d] =%18.7e" % (i, disp), file=stderr)
                print("   repl[%d] =%18.7e" % (i, repl), file=stderr)
                print("   epot[%d] =%18.7e" % (i, totl), file=stderr)
                print("sum_SR = %18.7e" % sum_SR, file=stderr)
                print("sum_14 = %18.7e" % sum_14, file=stderr)
                print("sum    = %18.7e" % sum_totl, file=stderr)

        # force 
        rij = conf.getDisplacements(self.ai, self.aj)
        invDist8 = dist ** (-8)
        f_i = 6 * (2*self.c12*invDist6 - self.c6) * np.transpose(rij) * invDist8
        f_j = -f_i
        for i in range(self.size):
            forceConf[self.ai[i],:] += np.transpose(f_i)[i,:]
            forceConf[self.aj[i],:] += np.transpose(f_j)[i,:]
        return energies

class coulombTerms (object):

    def __init__ (self, size):
        self.size = size
        self.ai = np.zeros(self.size, dtype=np.uint8)
        self.aj = np.zeros(self.size, dtype=np.uint8)
        self.qij = np.zeros(self.size)

    def setMember (self, i, ai, aj, qij):
        self.ai[i] = ai - 1
        self.aj[i] = aj - 1
        self.qij[i] = qij

    def calcForConf (self, conf, forceConf):
        if (self.size == 0):
            return 0
        dist = conf.getDistances (self.ai, self.aj)
        invDist1 = dist ** (-1)
        energies = 138.9354 * self.qij * invDist1
        # force 
        rij = conf.getDisplacements(self.ai, self.aj)
        f_i = (energies / (dist * dist)) * np.transpose(rij)
        f_j = -f_i
        for i in range(self.size):
            forceConf[self.ai[i],:] += np.transpose(f_i)[i,:]
            forceConf[self.aj[i],:] += np.transpose(f_j)[i,:]
        return energies

class MMCalculator (object):

    def __init__ (self):
        self.atomTerms = atomTerms(0)
        self.bondTerms = G96bondTerms(0)
        self.angleTerms = G96angleTerms(0)
        self.dihedralTerms = dihedralTerms(0)
        self.optDihedralTerms = optDihedralTerms(0)
        self.dihedralRestraintTerms = dihedralRestraintTerms(0)
        self.improperTerms = improperTerms(0)
        self.LJTerms = LJTerms(0)
        self.coulombTerms = coulombTerms(0)
        self.forceConf = None

    def createFromStpDictionary (self, stpDict):
        key   = 'atoms'
        ndofs = len(stpDict[key])
        natoms = ndofs
        self.atomTerms = atomTerms (natoms)
        for i in range(ndofs):
            info = stpDict[key][i]
            self.atomTerms.setMember (i, info['c6'], info['c12'], info['cs6'], info['cs12'], info['q'])
          
        key   = 'bonds'
        ndofs = len(stpDict[key][0])
        if (ndofs != 0):
            # get bond types and check if they are all the same
            bondtypes = [stpDict[key][0][i][2] for i in range(ndofs)]
            bondtype  = bondtypes[0]
            if bondtypes.count(bondtype) != ndofs:
                raise Exception("Not all bonds are of the same type. Check your .stp files.")
            if (bondtype == 1):
                self.bondTerms = harmonicBondTerms(ndofs)
            elif (bondtype == 2):
                self.bondTerms = G96bondTerms(ndofs)
            else:
                raise Exception("Bond type {} is not supported.".format(bondtype))
            for i in range(ndofs):
                self.bondTerms.setMember(i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][1][i], stpDict[key][2][i])

        key   = 'angles'
        ndofs = len(stpDict[key][0])
        if (ndofs != 0):
            # get angle types and check if they are all the same
            angletypes = [stpDict[key][0][i][3] for i in range(ndofs)]
            angletype  = angletypes[0]
            if angletypes.count(angletype) != ndofs:
                raise Exception("Not all angles are of the same type. Check your .stp files.")
            if (angletype == 1):
                self.angleTerms = harmonicAngleTerms(ndofs)
                for i in range(ndofs):
                    self.angleTerms.setMember (i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][1][i], stpDict[key][2][i] ) 
            elif (angletype == 2):
                self.angleTerms = G96angleTerms(ndofs) 
                for i in range(ndofs):
                    self.angleTerms.setMember (i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][1][i], stpDict[key][2][i] ) 
            elif (angletype == 5):
                self.angleTerms = RBAngleTerms(ndofs)
                for i in range(ndofs):
                    self.angleTerms.setMember (i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][1][i], stpDict[key][2][i], stpDict[key][3][i], stpDict[key][4][i] )
            else:
                raise Exception("Angle type {} is not supported.".format(angletype))

        key   = 'propers'
        ndofs = len(stpDict[key][0])
        if (ndofs != 0):
            # check if all types are the same
            dihtypes = [stpDict[key][0][i][4] for i in range(ndofs)]
            dihtype  = dihtypes[0]
            if dihtypes.count(dihtype) != ndofs:
                raise Exception("Not all proper dihedrals are of the same type. Check your .stp files.")
            if (dihtype == 1) or (dihtype == 9):
                # Periodic proper
                noptimized = len(stpDict['optdihedrals'])
                self.dihedralTerms    = dihedralTerms(ndofs - noptimized)
                self.optDihedralTerms = optDihedralTerms(noptimized)
                j = 0
                for i in range(ndofs):
                    if i not in stpDict['optdihedrals']:
                        self.dihedralTerms.setMember(j, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][0][i][3],
                                                     stpDict[key][1][i], stpDict[key][2][i], stpDict[key][3][i])
                        j += 1
                # now Opt
                for j,i in enumerate(stpDict['optdihedrals']):
                    self.optDihedralTerms.setMember(j, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][0][i][3])
                    
            elif (dihtype == 3):
                # Ryckaert-Bellemanns
                self.dihedralTerms = RyckaertBellemansDihedralTerms(ndofs)
                for i in range(ndofs):
                    self.dihedralTerms.setMember(i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][0][i][3],
                                                 stpDict[key][1][i], stpDict[key][2][i], stpDict[key][3][i], stpDict[key][4][i], stpDict[key][5][i], stpDict[key][6][i])
            elif (dihtype == 5):
                # Fourier
                self.dihedralTerms = FourierDihedralTerms(ndofs)
                for i in range(ndofs):
                    self.dihedralTerms.setMember(i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][0][i][3],
                                                 stpDict[key][1][i], stpDict[key][2][i], stpDict[key][3][i], stpDict[key][4][i])
                
        key   = 'impropers'
        ndofs = len(stpDict[key][0])
        if (ndofs != 0):
            # check if all types are the same
            dihtypes = [stpDict[key][0][i][4] for i in range(ndofs)]
            dihtype  = dihtypes[0]
            if dihtypes.count(dihtype) != ndofs:
                raise Exception("Not all improper dihedrals are of the same type. Check your .stp files.")
            if (dihtype == 2):
                self.improperTerms = improperTerms(ndofs)
                for i in range(ndofs):
                    self.improperTerms.setMember(i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][0][i][3], stpDict[key][1][i], stpDict[key][2][i] ) 
            elif (dihtype == 4):
                self.improperTerms = periodicImproperTerms(ndofs)
                for i in range(ndofs):
                    # note inverted order of k, phi
                    self.improperTerms.setMember(i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][0][i][3], stpDict[key][2][i], stpDict[key][1][i], stpDict[key][3][i] ) 
            else:
                raise Exception("Improper type {} is not supported.".format(dihtype))


        key   = 'restraints'
        ndofs = len(stpDict[key][0])
        self.dihedralRestraintTerms = dihedralRestraintTerms (ndofs)
        for i in range(ndofs):
            self.dihedralRestraintTerms.setMember(i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2], stpDict[key][0][i][3],
                stpDict[key][2][i], stpDict[key][1][i] ) 

        key   = 'nb'
        ndofs = len(stpDict[key][0])
        self.LJTerms = LJTerms (ndofs)
        self.coulombTerms = coulombTerms (ndofs)
        for i in range(ndofs):
            self.LJTerms.setMember(i, stpDict[key][0][i][0], stpDict[key][0][i][1], stpDict[key][0][i][2],
                stpDict[key][1][i], stpDict[key][2][i] ) 
            self.coulombTerms.setMember (i, stpDict[key][0][i][0], stpDict[key][0][i][1],
                stpDict[key][3][i] ) 

        # also initialize force configuration
        self.forceConf = np.zeros((natoms,3))

    def clearForces (self):
        self.forceConf = np.zeros(self.forceConf.shape)

    def forceNorm (self):
        return np.linalg.norm(self.forceConf)

    def getForceConf (self):
        return self.forceConf.copy()

    def setForces (self, forceconf):
        self.forceConf = forceconf.copy()

    def normalizeForces (self):
        norm = self.forceNorm()
        self.forceConf /= norm

    def applyForcesToConfWithFactor (self, conf, factor, force=None):
        if (force is None):
            conf += (factor * self.forceConf)
        else:
            conf += (factor * force)
            
    def setSingleDihedralRestraint (self, ai, aj, ak, al, phi_0, k):
        self.dihedralRestraintTerms =  dihedralRestraintTerms (1)
        self.dihedralRestraintTerms.setMember(0, ai, aj, ak, al, phi_0, k)

    def popDihedralRestraint (self):
        self.dihedralRestraintTerms.popMember()

    def pushDihedralRestraint (self, ai, aj, ak, al, phi_0, k):
        self.dihedralRestraintTerms.pushMember(ai, aj, ak, al, phi_0, k)

    def setLJParametersForAtom (self, i, cs6=None, cs12=None, mixtype='geometric'):
        if cs6 is not None:
            self.atomTerms.cs6[i] = cs6
        if cs12 is not None:
            self.atomTerms.cs12[i] = cs12
        # now update the LJ terms based on these new values
        self.LJTerms.setFromAtoms(self.atomTerms, mixtype)

    def setLJParametersForPair(self, i, cs6=None, cs12=None):
        self.LJTerms.setParameters(i, cs6, cs12)

    def setDihedralParameters (self, i, phi=None, k=None, m=None):
        self.dihedralTerms.setParameters(i,phi,k,m)

    def setOptDihedralParameters(self, m, phi=None, k=None):
        self.optDihedralTerms.setParameters(m, phi, k)

    def setDihedralParametersRyck (self, i, j, k):
        self.dihedralTerms.setParameters(i, j, k)

    # duplicates dihedral term "i" and puts copy at the end of the dihedral terms
    # returns index of the new dihedral
    def duplicateDihedral (self, i):
        return self.dihedralTerms.duplicateDihedral(i)

    def calcForConf (self, conf, removeRestraintsFromTotal=False, debug=False):
        self.clearForces()
        outputDict = {}
        outputDict['total'] = 0.0
        outputDict['bonds'] = np.sum(self.bondTerms.calcForConf(conf, self.forceConf))
        outputDict['angles'] = np.sum(self.angleTerms.calcForConf(conf, self.forceConf))
        outputDict['propers'] = np.sum(self.dihedralTerms.calcForConf(conf, self.forceConf)) + np.sum(self.optDihedralTerms.calcForConf(conf, self.forceConf))
        outputDict['impropers'] = np.sum(self.improperTerms.calcForConf(conf, self.forceConf))
        outputDict['lj'] = np.sum(self.LJTerms.calcForConf(conf, self.forceConf, debug=debug))
        outputDict['coulomb'] = np.sum(self.coulombTerms.calcForConf(conf, self.forceConf))
        outputDict['restraints'] = np.sum(self.dihedralRestraintTerms.calcForConf(conf, self.forceConf))
        outputDict['total'] = sum(outputDict.values())
        if (removeRestraintsFromTotal):
            outputDict['total'] -= outputDict['restraints']
        return (outputDict, self.forceConf.copy())

    def calcForEnsemble (self, ens, shiftToZero=True, removeRestraintsFromTotal=True, debug=False):
        outputDict = {}
        outputDict['bonds'] = []
        outputDict['angles'] = []
        outputDict['propers'] = []
        outputDict['restraints'] = []
        outputDict['impropers'] = []
        outputDict['lj'] = []
        outputDict['coulomb'] = []
        outputDict['total'] = []
        for conf in ens:
            thisMember = self.calcForConf (conf, removeRestraintsFromTotal, debug=debug)[0]
            for key in outputDict:
                outputDict[key].append(thisMember[key])
        if (shiftToZero):
            minim = min(outputDict['total'])
            for key in outputDict:
                for i in range(len(outputDict[key])):
                    outputDict[key][i] -= minim
        return outputDict

    def calcForEnsembleAndSaveToFile (self, ens, fn, shiftToZero=True, removeRestraints=True, saveOnlyTotal=True, debug=False):
        data = self.calcForEnsemble (ens, shiftToZero, removeRestraints, debug)
        if (saveOnlyTotal):
            np.savetxt(fn, data['total'])
        else:
            fp = open(fn, "w")
            if (shiftToZero):
                fp.write("# Note that all potentials are shifted so as to keep the total potential above 0.\n")
            fp.write("@ xaxis label \"Configuration\"\n")
            fp.write("@ yaxis label \"Energy (kJ/mol)\"\n")
            # begin with total
            fp.write("@ s0 legend \"total\"\n")
            i = 1
            for key in data:
                if key != 'total' and not (saveOnlyTotal):
                    fp.write("@ s%d legend \"%s\"\n" % (i,key))
                    i += 1
            for i in range(ens.size()):
                fp.write("%10d" % (i+1))
                fp.write("%18.7e" % data['total'][i])
                for key in data:
                    if key != 'total' and not (saveOnlyTotal):
                        fp.write("%18.7e" % data[key][i])
                fp.write("\n")
            fp.close()
