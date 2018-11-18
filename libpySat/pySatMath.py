import numpy as np

def RotationMatrix3DX( theta: float):
    """
    system rotation matrix around X axis
    :param theta: angle in radians
    :return:
    """
    ret = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)

    ret[0, 0] = 1.0
    ret[0, 1] = 0.0
    ret[0, 2] = 0.0
    ret[1, 0] = 0.0;
    ret[1, 1] = c
    ret[1, 2] = s
    ret[2, 0] = 0.0
    ret[2, 1] = -s
    ret[2, 2] = c
    return ret

def RotationMatrix3DY( theta: float):
    """
    system rotation matrix around Y axis
    :param theta: angle in radians
    :return:
    """
    ret = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)

    ret[0, 0] = c
    ret[0, 1] = 0.0
    ret[0, 2] = -s
    ret[1, 0] = 0.0;
    ret[1, 1] = 1.0
    ret[1, 2] = 0.0
    ret[2, 0] = s
    ret[2, 1] = 0.0
    ret[2, 2] = c
    return ret

def RotationMatrix3DZ( theta: float):
    """
    system rotation matrix around Z axis
    :param theta: angle in radians
    :return:
    """
    ret = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)

    ret[0, 0] = c
    ret[0, 1] = s
    ret[0, 2] = 0.0
    ret[1, 0] = -s
    ret[1, 1] = c
    ret[1, 2] = 0.0
    ret[2, 0] = 0.0
    ret[2, 1] = 0.0
    ret[2, 2] = 1.0
    return ret

def InterpolateLagrange(epochWanted: float, nodes: np.array, values: np.array):
    """
    Lagrange interpolation
    :param epochWanted:
    :param nodes:
    :param values:
    :return:
    """
    n= len(nodes)
    ret=0.0
    for i in range(0,n):
        tmp =1.0
        for j in range(0,n):
            if (i==j):
                tmp*=1
            else:
                tmp*= ( epochWanted - nodes[j] ) / ( nodes[i] - nodes[j] )
        ret += tmp* values[i]
    return ret

#def InterpolateLagrangeD(epochWanted: float, nodes: np.array, values: np.array):
#    """
#    Lagrange interpolation returning the time derivative.
#    :param epochWanted:
#    :param nodes:
#    :param values:
#    :return:
#    """
#
#template<class T, class X> X
#InterpolateD ( ublas::vector<T> timenodes, ublas::vector<X> y, T timewanted )
#{
#    int nx = timenodes.size();
#    X ret = 0.0;
#
#    for ( int i = 0; i < nx; i++ )
#    {
#        double tmp = 0.0;
#        for ( int j = 0; j < nx; j++ )
#        {
#            if ( i == j ) continue;
#            double tmp2 = 1.0;
#            for ( int jj = 0; jj < nx; jj++ )
#            {
#                if ( i == jj ) continue;
#                if ( j != jj ) tmp2 *= ( timewanted - timenodes[jj] ) / ( timenodes[i] - timenodes[jj] );
#                if ( j == jj ) tmp2 *= 1.0 / ( timenodes[i] - timenodes[j] );
#            }
#            tmp += tmp2;
#        }
#        ret += tmp * y ( i );
#    }
#    return ret;
#}