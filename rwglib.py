from scipy.linalg import solve
import numpy as np
import numpy.matlib

##IMPMET Standard impedance matrix (metal surface)
##Returns the complex impedance matrix
##           [EdgesTotal x EdgesTotal]
##Uses 9 integration points for every triangle (barycentric subdivision)
##The impedance matrix is calculated as a sum of the contributions due to separate triangles (similar to the "face-pair" method). See Appendix B for a detailed algorithm.
##A 9-point quadrature is used for all integrals, including the self-coupling terms. The alternative source code with the analytical approximation of the self-coupling terms is given in Appendix B. The difference between two methods is not significant. 

def impmet(EdgesTotal,TrianglesTotal,EdgeLength,K,Center,Center_,TrianglePlus,TriangleMinus,RHO_P,RHO_M,RHO__Plus,RHO__Minus,FactorA,FactorFi):
    # Memory allocation
    Z   = np.zeros ((EdgesTotal,EdgesTotal))+1j*np.zeros((EdgesTotal,EdgesTotal)) 
    print("Impedance matrice's size: ",Z.shape)
    # Loop over integration triangles
    for p_ in range(TrianglesTotal):
        Plus  = [i for i,j in enumerate(np.array(TrianglePlus)-p_) if j == 0]
        Minus  = [i for i,j in enumerate(np.array(TriangleMinus)-p_) if j == 0]

        D = Center_ - np.transpose(np.tile([Center[:,p_]],(TrianglesTotal,9,1))) #[3 9 TrianglesTotal]     

        R=np.sqrt(sum(D*D))[np.newaxis]                              #[1 9 TrianglesTotal]
        g=np.exp(-K*R)/R                                #[1 9 TrianglesTotal]

        gP=g[:,:,TrianglePlus]                         #[1 9 EdgesTotal]
        gM=g[:,:,TriangleMinus]                        #[1 9 EdgesTotal]

        Fi= (np.sum(gP,axis = 1) - np.sum(gM,axis = 1))#[1 1 EdgesTotal]
        ZF= FactorFi*Fi.T         #[EdgesTotal 1]

        for n in Plus:
            #n = Plus[k]
            RP = np.tile([RHO__Plus[:,:,n]],(EdgesTotal,1,1)).transpose(1,2,0)  #[3 9 EdgesTotal]
            A=(np.sum(gP*np.sum(RP*RHO_P,axis=0),axis=1)+np.sum(gM*np.sum(RP*RHO_M,axis=0),axis=1))
            Z1= FactorA*A.T
            Z[:,n]=Z[:,n] + (EdgeLength[n]*(Z1+ZF)).reshape(EdgesTotal)

        for n in Minus:
            #n = Minus[k]
            RP = np.tile([RHO__Minus[:,:,n]],(EdgesTotal,1,1)).transpose(1,2,0)  #[3 9 EdgesTotal]
            A=(np.sum(gP*np.sum(RP*RHO_P,axis=0),axis=1)+np.sum(gM*np.sum(RP*RHO_M,axis=0),axis=1))
            Z1= FactorA*A.T
            Z[:,n]=Z[:,n] + (EdgeLength[n]*(Z1-ZF)).reshape(EdgesTotal)
    return Z
