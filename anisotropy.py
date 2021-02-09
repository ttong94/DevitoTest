import numpy as np
import math


# __all__ =  # look up what's this mean?

# C = np.array([[C11, C12, C13, C14, C15, C16],
#               [C12, C22, C23, C24, C25, C26],
#               [C13, C23, C33, C34, C35, C36],
#               [C14, C24, C34, C44, C45, C46],
#               [C15, C25, C35, C45, C55, C56],
#               [C16, C26, C36, C46, C56, C66]])
"""
S = np.array([[A, B, C, D, E, F],
              [B, G, H, I, J, K],
              [C, H, L, M, N, O],
              [D, I, M, P, Q, R],
              [E, J, N, Q, S, T],
              [F, K, O, R, T, U]])

so C
"""

class StiffnessMatrix:
    """
    by default it is an isotropic

    """

    def __init__(self, C11 = 20.37, C12 = 12.3, C13 = 12.3, C14 = 0, C15 = 0, C16 = 0,
                                   C22 = 20.37, C23 = 12.3, C24 = 0, C25 = 0, C26 = 0,
                                                C33 = 20.37, C34 = 0, C35 = 0, C36 = 0,
                                                            C44 = 4.035, C45 = 0, C46 = 0,
                                                                        C55 = 4.035, C56 = 0,
                                                                                      C66 = 4.035):
        self.data = np.array([ [C11, C12, C13, C14, C15, C16],
                          [C12, C22, C23, C24, C25, C26],
                          [C13, C23, C33, C34, C35, C36],
                          [C14, C24, C34, C44, C45, C46],
                          [C15, C25, C35, C45, C55, C56],
                          [C16, C26, C36, C46, C56, C66]])


    def rotate(self, Q):
        """
        Q is the rotation matrix, can be obtained from rotator.getMatrix.
        reference: http://www.continuummechanics.org/coordxforms.html
        """
#         return np.matmul(np.matmul(np.matmul(np.matmul(Q.T,Q.T), self.get4DTensor), Q), Q)
        return np.matmul( np.matmul(np.matmul(np.matmul(Q.T,Q.T), self.get4DTensor()), Q), Q)

    
    def get4DTensor(self, compl=False):
        """Convert from Voigt to full tensor notation 
           Convert from the 6*6 elastic constants matrix to 
           the 3*3*3*3 tensor representation. Recoded from 
           the Fortran implementation in DRex. Use the optional 
           argument "compl" for the elastic compliance (not 
           stiffness) tensor to deal with the multiplication 
           of elements needed to keep the Voigt and full 
           notation consistant.
        """
#         print(self.data)
        from sympy import Array
        cij_mat = self.data
        cij_tens = np.zeros((3,3,3,3))
        m2t = np.array([[0,5,4],[5,1,3],[4,3,2]])
        if compl:
            cij_mat = cij_mat / np.array([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                          [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                          [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                          [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                          [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                          [2.0, 2.0, 2.0, 4.0, 4.0, 4.0]])
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        cij_tens[i,j,k,l] = cij_mat[m2t[i,j],m2t[k,l]]

        return cij_tens


class StiffnessMatrixVTI(StiffnessMatrix):
    def __init__(self, data):
        self.data = data
        self.epsilon = self.getEpsilon()
        self.gamma = self.getGamma()
        self.delta = self.getDelta()


    def getEpsilon(self):
        return (self.data[0,0] - self.data[2,2]) / (2 * self.data[2,2])


    def getGamma(self):
        return (self.data[5,5] - self.data[3,3]) / (2 * self.data[3,3])


    def getDelta(self):
        return (((self.data[0,2] + self.data[3,3])**2 - (self.data[2,2] - self.data[3,3])**2)
                /(2 * self.data[2,2] * (self.data[2,2] - self.data[3,3])))




class Rotator:
    """
    reference:
    Page 9.
    https://ocw.mit.edu/courses/materials-science-and-engineering/3-11-mechanics-of-materials-fall-1999/modules/MIT3_11F99_trans.pdf
    """

    def __init__(self, phi = 0, psi = 0, theta = 0):
        """

        """
        self.phi = phi
        self.psi = psi
        self.theta = theta

    def getMatrix(self):
        phi_matrix  = np.array([[math.cos(self.phi), math.sin(self.phi), 0],
                               [-math.sin(self.phi), math.cos(self.phi), 0],
                               [0                  , 0                 , 1]])
        psi_matrix  = np.array([[math.cos(self.psi), math.sin(self.psi), 0],
                               [-math.sin(self.psi), math.cos(self.psi), 0],
                               [0                  , 0                 , 1]])
        theta_matrix = np.array([[1, 0                   ,  0                   ],
                                 [0, math.cos(self.theta),  math.sin(self.theta)],
                                 [0, -math.sin(self.theta), math.cos(self.theta)]])

        return np.matmul(np.matmul(psi_matrix, theta_matrix), phi_matrix)














