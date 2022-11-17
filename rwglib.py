from scipy.linalg import solve
import numpy as np
import numpy.matlib

##IMPMET Standard impedance matrix (metal surface)
##Returns the complex impedance matrix
##           [EdgesTotal x EdgesTotal]
##Uses 9 integration points for every triangle (barycentric subdivision)
##The impedance matrix is calculated as a sum of the contributions due to separate triangles (similar to the "face-pair" method). See Appendix B for a detailed algorithm.
##A 9-point quadrature is used for all integrals, including the self-coupling terms. The alternative source code with the analytical approximation of the self-coupling terms is given in Appendix B. The difference between two methods is not significant. 

def impmet(rwgmesh,K,moment,FactorA,FactorFi):
    # Memory allocation
    Z   = np.zeros ((rwgmesh.EdgesTotal,rwgmesh.EdgesTotal))+1j*np.zeros((rwgmesh.EdgesTotal,rwgmesh.EdgesTotal)) 
    print("Impedance matrice's size: ",Z.shape)
    # Loop over integration triangles
    for p_ in range(rwgmesh.TrianglesTotal):
        Plus  = [i for i,j in enumerate(np.array(rwgmesh.TrianglePlus)-p_) if j == 0]
        Minus  = [i for i,j in enumerate(np.array(rwgmesh.TriangleMinus)-p_) if j == 0]

        D = moment.Center_ - np.transpose(np.tile([rwgmesh.Center[:,p_]],(rwgmesh.TrianglesTotal,9,1))) #[3 9 TrianglesTotal]     

        R=np.sqrt(sum(D*D))[np.newaxis]                              #[1 9 TrianglesTotal]
        g=np.exp(-K*R)/R                                #[1 9 TrianglesTotal]

        gP=g[:,:,rwgmesh.TrianglePlus]                         #[1 9 EdgesTotal]
        gM=g[:,:,rwgmesh.TriangleMinus]                        #[1 9 EdgesTotal]

        Fi= (np.sum(gP,axis = 1) - np.sum(gM,axis = 1))#[1 1 EdgesTotal]
        ZF= FactorFi*Fi.T         #[EdgesTotal 1]

        for n in Plus:
            #n = Plus[k]
            RP = np.tile([moment.RHO__Plus[:,:,n]],(rwgmesh.EdgesTotal,1,1)).transpose(1,2,0)  #[3 9 EdgesTotal]
            A=(np.sum(gP*np.sum(RP*moment.RHO_P,axis=0),axis=1)+np.sum(gM*np.sum(RP*moment.RHO_M,axis=0),axis=1))
            Z1= FactorA*A.T
            Z[:,n]=Z[:,n] + (rwgmesh.EdgeLength[n]*(Z1+ZF)).reshape(rwgmesh.EdgesTotal)

        for n in Minus:
            #n = Minus[k]
            RP = np.tile([moment.RHO__Minus[:,:,n]],(rwgmesh.EdgesTotal,1,1)).transpose(1,2,0)  #[3 9 EdgesTotal]
            A=(np.sum(gP*np.sum(RP*moment.RHO_P,axis=0),axis=1)+np.sum(gM*np.sum(RP*moment.RHO_M,axis=0),axis=1))
            Z1= FactorA*A.T
            Z[:,n]=Z[:,n] + (rwgmesh.EdgeLength[n]*(Z1-ZF)).reshape(rwgmesh.EdgesTotal)
    return Z

##RWG1: Geometry calculations - all Chapters
##Uses the structure mesh file, e.g. platefine.mat, as an input.
##Creates the RWG edge element for every inner edge of the structure. The total number of elements is EdgesTotal.
##Outputs the following arrays:
##Edge first node number
##        Edge_(1,1:EdgesTotal)
##Edge second node number
##        Edge_(2,1:EdgesTotal)
##Plus triangle number
##        TrianglePlus(1:EdgesTotal)
##Minus triangle number
##        TriangleMinus(1:EdgesTotal)
##Edge length
##        EdgeLength(1:EdgesTotal)
##Edge element indicator
##        EdgeIndicator(1:EdgesTotal)
##Also outputs areas and midpoints of separate triangles:
##Triangle area
##        Area(1:TrianglesTotal)
##Triangle center
##        Center(1:TrianglesTotal)      
##This script may handle surfaces with T-junctions including monopoles over various metal surfaces and certain metal meshes

class RWGmesh:
    def __init__(self, p,t):
        
        [s1,s2] = p.shape
        if s1 == 2:
            print("Add third dimension")
            p = np.append(p,np.zeros((1,s2)),axis=0) 
        
        self.TrianglesTotal=t.shape[1]
        print("Total number of triangles is ",self.TrianglesTotal)
        
        # Find areas of separate triangles
        self.Area = np.zeros(self.TrianglesTotal)
        self.Center = np.zeros((3,self.TrianglesTotal))
        for m in range(self.TrianglesTotal):
            N = t[0:3,m] - 1
            #print(N)
            Vec1 = p[:,N[0]] - p[:,N[1]]
            Vec2 = p[:,N[2]] - p[:,N[1]]
            self.Area[m] = np.linalg.norm(np.cross(Vec1,Vec2))/2
            self.Center[:,m] = 1/3*sum(np.transpose(np.array(p[:,N])))
        
        # Find all edge elements "Edge_" with at least two adjacent triangles
        n = 0
        Edge_l = []
        self.TrianglePlus = []
        self.TriangleMinus = []
        for m in range(self.TrianglesTotal):
            N = t[0:3,m] - 1 
            for k in range(m+1,self.TrianglesTotal):
                M = t[0:3,k] - 1      
                a = 1 - np.array([1 if i else 0 for i in [all(N-M[0]),all(N-M[1]),all(N-M[2])]])
                if(sum(a)==2): #triangles m and k have common edge
                    Edge_l.append([M[i] for i,j in enumerate(a) if j]); 
                    self.TrianglePlus.append(m);
                    self.TriangleMinus.append(k); 
                    n = n + 1;
        self.Edge_ = np.transpose(np.array(Edge_l))
        self.Edge__ = np.transpose(np.array([[i[1], i[0]] for i in Edge_l]))
        self.EdgesTotal = len(Edge_l)
        print("Total number of edges is ",self.EdgesTotal)
        
        # find edge length
        self.EdgeLength = np.zeros(self.EdgesTotal)
        for m in range(self.EdgesTotal):
            self.EdgeLength[m]=np.linalg.norm(p[:,self.Edge_[0,m]]-p[:,self.Edge_[1,m]])

##RWG2: Geometry calculations - all Chapters
##Uses the mesh file from RWG1, mesh1.mat, as an input.
##Creates the following parameters of the RWG edge elements:
##Position vector rho_c_plus from the free vertex of the "plus" triangle to its center
##                               RHO_Plus(1:3,1:EdgesTotal)
##Position vector rho_c_minus from the center of the "minus" triangle to its free vertex
##                               RHO_Minus(1:3,1:EdgesTotal)
##In addition to these parameters creates the following arrays for nine subtriangles (barycentric subdivision):
##Midpoints of nine subtriangles
##                               Center_(1:3,1:9,1:TrianglesTotal)  
##Position vectors rho_c_plus from the free vertex of the "plus" triangle to nine subtriangle midpoints
##                               RHO__Plus(1:3,1:9,1:EdgesTotal)
##Position vectors rho_c_minus from nine subtriangle midpoints to the free vertex of the "minus" triangle
##                               RHO__Minus(1:3,1:9,1:EdgesTotal)
##See Rao, Wilton, Glisson, IEEE Trans. Antennas and Propagation, vol. AP-30, No 3, pp. 409-418, 1982.

class RWGmesh2:
    def __init__(self, mesh, t, p):
        IMT=[]
        for m in range(mesh.TrianglesTotal):
            n1 = t[0,m] - 1
            n2 = t[1,m] - 1
            n3 = t[2,m] - 1
            M = mesh.Center[:,m]
            r1 = p[:,n1]
            r2 = p[:,n2]
            r3 = p[:,n3]
            r12 = r2-r1
            r23 = r3-r2
            r13 = r3-r1
            C1 = r1 + (1/3)*r12
            C2 = r1+(2/3)*r12
            C3 = r2+(1/3)*r23
            C4 = r2+(2/3)*r23
            C5 = r1+(1/3)*r13
            C6 = r1+(2/3)*r13
            a1 = np.transpose([1/3*(C1+C5+r1)])
            a2 = np.transpose([1/3*(C1+C2+M)])
            a3 = np.transpose([1/3*(C2+C3+r2)])
            a4 = np.transpose([1/3*(C2+C3+M)])
            a5 = np.transpose([1/3*(C3+C4+M)])
            a6 = np.transpose([1/3*(C1+C5+M)])
            a7 = np.transpose([1/3*(C5+C6+M)])
            a8 = np.transpose([1/3*(C4+C6+M)])
            a9 = np.transpose([1/3*(C4+C6+r3)])
            if m == 0:
                self.Center_ = [np.concatenate((a1,a2,a3,a4,a5,a6,a7,a8,a9), axis=1)]
            else:
                self.Center_ = np.concatenate((self.Center_,[np.concatenate((a1,a2,a3,a4,a5,a6,a7,a8,a9), axis=1)]), axis = 0)

        self.Center_ = self.Center_.transpose(1,2,0)
        
        # PLUS
        for m in range(mesh.EdgesTotal):
            NoPlus=mesh.TrianglePlus[m]
            n1 = t[0,NoPlus] - 1
            n2 = t[1,NoPlus] - 1
            n3 = t[2,NoPlus] - 1
            if((n1 != mesh.Edge_[0,m])&(n1 != mesh.Edge_[1,m])): 
                NODE=n1
            if((n2 != mesh.Edge_[0,m])&(n2 != mesh.Edge_[1,m])): 
                NODE=n2
            if((n3 != mesh.Edge_[0,m])&(n3 != mesh.Edge_[1,m])): 
                NODE=n3
            FreeVertex=p[:,NODE]

            if m == 0:
                self.RHO_Plus   = [mesh.Center[:,NoPlus]-FreeVertex]
                #Nine rho's of the "plus" triangle
                self.RHO__Plus  = [self.Center_[:,:,NoPlus] - np.transpose(np.tile([FreeVertex],(9,1)))]
            else:
                self.RHO_Plus   = np.concatenate((self.RHO_Plus,[mesh.Center[:,NoPlus]-FreeVertex]), axis = 0)
                self.RHO__Plus  = np.concatenate((self.RHO__Plus,[self.Center_[:,:,NoPlus] - np.transpose(np.tile([FreeVertex],(9,1)))]), axis = 0)

        self.RHO_Plus = self.RHO_Plus.transpose(1,0)
        self.RHO__Plus = self.RHO__Plus.transpose(1,2,0)
        
        m = 0
        NoPlus=mesh.TrianglePlus[m]
        n1 = t[0,NoPlus] - 1
        n2 = t[1,NoPlus] - 1
        n3 = t[2,NoPlus] - 1
        if((n1 != mesh.Edge_[0,m])&(n1 != mesh.Edge_[1,m])): 
            NODE=n1
        if((n2 != mesh.Edge_[0,m])&(n2 != mesh.Edge_[1,m])): 
            NODE=n2
        if((n3 != mesh.Edge_[0,m])&(n3 != mesh.Edge_[1,m])): 
            NODE=n3
        FreeVertex=p[:,NODE]
        
        # MINUS
        for m in range(mesh.EdgesTotal):
            NoMinus=mesh.TriangleMinus[m]
            n1 = t[0,NoMinus] - 1
            n2 = t[1,NoMinus] - 1
            n3 = t[2,NoMinus] - 1
            if((n1 != mesh.Edge_[0,m])&(n1 != mesh.Edge_[1,m])): 
                NODE=n1
            if((n2 != mesh.Edge_[0,m])&(n2 != mesh.Edge_[1,m])): 
                NODE=n2
            if((n3 != mesh.Edge_[0,m])&(n3 != mesh.Edge_[1,m])): 
                NODE=n3
            FreeVertex=p[:,NODE]

            if m == 0:
                self.RHO_Minus   = [-mesh.Center[:,NoMinus] + FreeVertex]
                #Nine rho's of the "plus" triangle
                self.RHO__Minus  = [-self.Center_[:,:,NoMinus] + np.transpose(np.tile([FreeVertex],(9,1)))]
            else:
                self.RHO_Minus   = np.concatenate((self.RHO_Minus,[-mesh.Center[:,NoMinus] + FreeVertex]), axis = 0)
                self.RHO__Minus  = np.concatenate((self.RHO__Minus,[-self.Center_[:,:,NoMinus] + np.transpose(np.tile([FreeVertex],(9,1)))]), axis = 0)

        self.RHO_Minus = self.RHO_Minus.transpose(1,0)
        self.RHO__Minus = self.RHO__Minus.transpose(1,2,0)
        
        for m in range(mesh.EdgesTotal):
            if m == 0:
                self.RHO_P = [np.transpose(np.tile([self.RHO_Plus[:,m]],(9,1)))]   #[3 9 EdgesTotal]
                self.RHO_M = [np.transpose(np.tile([self.RHO_Minus[:,m]],(9,1)))]  #[3 9 EdgesTotal]
            else:
                self.RHO_P = np.concatenate((self.RHO_P,[np.transpose(np.tile([self.RHO_Plus[:,m]],(9,1)))]), axis = 0)   #[3 9 EdgesTotal]
                self.RHO_M = np.concatenate((self.RHO_M,[np.transpose(np.tile([self.RHO_Minus[:,m]],(9,1)))]), axis = 0)  #[3 9 EdgesTotal]
                
        self.RHO_P = self.RHO_P.transpose(1,2,0)
        self.RHO_M = self.RHO_M.transpose(1,2,0)
