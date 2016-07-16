import numpy as np
from warp import top

def FRBT(beta=5.0, alpha=0.0):
    """ 
        Transforms a matched flat beam to a round 'magnetized' beam.
    """

    gamma = (1. - alpha**2) / beta

    R = np.zeros([6,6],dtype='float64')
    R[0,0] = 1. + alpha
    R[0,1] = beta
    R[0,2] = 1. - alpha
    R[0,3] = -beta

    R[1,0] = -gamma
    R[1,1] = 1. - alpha
    R[1,2] = gamma
    R[1,3] = 1. + alpha

    R[2,0] = 1. - alpha
    R[2,1] = -beta
    R[2,2] = 1. + alpha
    R[2,3] = beta

    R[3,0] = gamma
    R[3,1] = 1. + alpha
    R[3,2] = -gamma
    R[3,3] = 1. - alpha

    R[4,4] = 2.
    R[5,5] = 2.

    R = 0.5 * R
    
    x = {}

    norm = {}
    for i in range(6):
        for j in range(6):
            norm[i,j] = 1.0
    norm[0,1] = norm[0,3] = norm[2,1] = norm[2,3] = 1./top.pgroup.uzp
    norm[1,0] = norm[1,2] = top.pgroup.uzp
    norm[3,0] = norm[3,2] = top.pgroup.uzp

    x = {}
    x[0] = np.copy(top.pgroup.xp)
    x[1] = np.copy(top.pgroup.uxp)
    x[2] = np.copy(top.pgroup.yp)
    x[3] = np.copy(top.pgroup.uyp)
    x[4] = np.copy(top.pgroup.zp)
    x[5] = np.copy(top.pgroup.uzp)

    print x[0].shape

    holding = []

    for i in range(6):
        val = 0
        for j in range(6):
            val += R[i,j] * x[j] * norm[i,j]
        
        holding.append(val)


    top.pgroup.xp = holding[0]
    top.pgroup.uxp = holding[1]
    top.pgroup.yp = holding[2]
    top.pgroup.uyp = holding[3]
    top.pgroup.zp = holding[4]
    top.pgroup.uzp = holding[5]
    
#     print "Transform!"


def IFRBT(beta=5.0, alpha=0.0):
    """ 
        Transforms a round 'magnetized' beam back to a flat beam.
    """

    gamma = (1. - alpha**2) / beta

    R = np.zeros([6,6],dtype='float64')
    R[0,0] = 1. + alpha
    R[0,1] = beta
    R[0,2] = alpha - 1.
    R[0,3] = beta

    R[1,0] = -gamma
    R[1,1] = 1. - alpha
    R[1,2] = -gamma
    R[1,3] = -1. - alpha

    R[2,0] = alpha - 1.
    R[2,1] = beta
    R[2,2] = 1. + alpha
    R[2,3] = beta

    R[3,0] = -gamma
    R[3,1] = -1. - alpha
    R[3,2] = -gamma
    R[3,3] = 1. - alpha

    R[4,4] = 2.
    R[5,5] = 2.

    R = 0.5 * R
    
    x = {}

    norm = {}
    for i in range(6):
        for j in range(6):
            norm[i,j] = 1.0
    norm[0,1] = norm[0,3] = norm[2,1] = norm[2,3] = 1./top.pgroup.uzp
    norm[1,0] = norm[1,2] = top.pgroup.uzp
    norm[3,0] = norm[3,2] = top.pgroup.uzp

    x = {}
    x[0] = np.copy(top.pgroup.xp)
    x[1] = np.copy(top.pgroup.uxp)
    x[2] = np.copy(top.pgroup.yp)
    x[3] = np.copy(top.pgroup.uyp)
    x[4] = np.copy(top.pgroup.zp)
    x[5] = np.copy(top.pgroup.uzp)

    print x[0].shape

    holding = []

    for i in range(6):
        val = 0
        for j in range(6):
            val += R[i,j] * x[j] * norm[i,j]
        
        holding.append(val)


    top.pgroup.xp = holding[0]
    top.pgroup.uxp = holding[1]
    top.pgroup.yp = holding[2]
    top.pgroup.uyp = holding[3]
    top.pgroup.zp = holding[4]
    top.pgroup.uzp = holding[5]

#     print "Inverse Transform!"