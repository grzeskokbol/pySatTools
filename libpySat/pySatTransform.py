

import datetime
import numpy as np
import libpySat as pySat
from astropy import _erfa as erfa
import warnings
from scipy import interpolate


class Transform:
    def __init__(self):

        self.epochSave=datetime.datetime.now()
        self.rotSave =np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]),dtype=float)
        self.rotSaveDot=self.rotSave
        self.iers_tab = pySat.GetIersTable()
        self.ut1utc=0.0
        self.ut1 = 0.0
        self.fut1utc, self.fxp, self.fyp, self.fdX, self.fdY = self.__getUtPoNut()
        self.precNut=pySat.TransformPrecessionNutation(self.fdX,self.fdY)
        self.polMot = pySat.TransformPolarMotion(self.fxp,self.fyp)
        self.era=pySat.TransformERA(self.fut1utc)

    def __getUtPoNut(self):
        """
        Gets UT1-UTC, Polar Motion and dX, dY from IERS A or IERS B and
        assigns to the class' variables
        :return:
        """

        fut1utc = interpolate.interp1d( self.iers_tab['MJD'], self.iers_tab['UT1_UTC'])
        fxp = interpolate.interp1d(self.iers_tab['MJD'], self.iers_tab['PM_x'])
        fyp = interpolate.interp1d( self.iers_tab['MJD'], self.iers_tab['PM_y'])
        fdX = interpolate.interp1d( self.iers_tab['MJD'], self.iers_tab['dX_2000A'])
        fdY = interpolate.interp1d( self.iers_tab['MJD'], self.iers_tab['dY_2000A'])


        return fut1utc, fxp, fyp, fdX,fdY

    def __RotationMatrixTRF2CRF(self,epoch:datetime.datetime):
        """
        Computes rotation matrix from ITRF to ICRF
        :param epoch:
        :return: 3x3 rotation matrix = Q R W
        """

        if (epoch != self.epochSave):
            ret=np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]),dtype=float)
            ret =self.precNut.getMatrix_PrecessionNutation(epoch)
            #print('PrecNutation',ret)
            ret =np.matmul(ret,self.era.getMatrix_Era(epoch))
            #print('ERAmatrix',getMatrix_Era(epoch))
            ret =np.matmul(ret,self.polMot.getMatrix_PolarMotion(epoch))
            #print('PolarMotion',getMatrix_PolarMotion(epoch))
            epochSave = epoch

            self.rotSave = ret
            self.epochSave=epoch
            #print("RotationMatrix ",ret)
            return self.rotSave
        else:
            return self.rotSave

    def __RotationMatrixDotTRF2CRF(self, epoch: datetime.datetime):
        """
        Computes rotation matrix from ITRF to ICRF
        :param epoch:
        :return: 3x3 rotation matrix = Q R W
        """

        rq = self.precNut.getMatrix_PrecessionNutation(epoch)
        rqdot = self.precNut.getMatrix_PrecessionNutationDot(epoch)
        # print('PrecNutation',ret)
        rr =  self.era.getMatrix_EraDot(epoch)
        rrdot =  self.era.getMatrix_EraDot(epoch)
        # print('ERAmatrix',getMatrix_Era(epoch))
        rw=self.polMot.getMatrix_PolarMotion(epoch)
        rwdot = self.polMot.getMatrix_PolarMotionDot(epoch)
        # print('PolarMotion',getMatrix_PolarMotion(epoch))
        r1 =np.matmul(rqdot,rr)
        r2 = np.matmul(rq, rrdot)
        r3 = np.matmul(rq, rr)

        self.rotSaveDot=np.matmul(r1,rw)+np.matmul(r2,rw)+np.matmul(r3,rwdot)
        return self.rotSaveDot

    def __RotationMatrixCRF2TRF(self,epoch:datetime.datetime):
        """
        Computes rotation matrix from ICRF to ITRF
        :param epoch:
        :return:
        """

        return self.__RotationMatrixTRF2CRF(epoch).transpose()

    def __RotationMatrixDotCRF2TRF(self,epoch:datetime.datetime):
        """
        Computes rotation matrix from ICRF to ITRF
        :param epoch:
        :return:
        """

        return self.__RotationMatrixDotTRF2CRF(epoch).transpose()

    def TRF2CRF(self,epoch:datetime.datetime,stateXyz: np.array):
        """

        :param epoch:
        :param stateXyz:
        :return:
        """
        posXyz=np.array([0 , 0, 0],dtype=float)
        velXyz = np.array([0,0 ,0], dtype=float)

        if (len(stateXyz)==3) :
            posXyz=stateXyz

        else:
            posXyz=stateXyz[0:3]
            velXyz=stateXyz[3:6]

        ret = np.array([0, 0, 0, 0, 0, 0], dtype=float)
        ret[0:3]=self.__RotationMatrixTRF2CRF(epoch).dot(posXyz)
        ret[3:6]=self.__RotationMatrixTRF2CRF(epoch).dot(velXyz) +self.__RotationMatrixDotTRF2CRF(epoch).dot(posXyz)

        return ret

    def CRF2TRF(self,epoch:datetime.datetime,stateXyz: np.array):
        """

        :param epoch:
        :param posXyz:
        :return:
        """
        posXyz = np.array([0, 0, 0], dtype=float)
        velXyz = np.array([0, 0, 0], dtype=float)

        if (len(stateXyz) == 3):
            posXyz = stateXyz

        else:
            posXyz = stateXyz[0:3]
            velXyz = stateXyz[3:6]

        ret = np.array([0, 0, 0, 0, 0, 0], dtype=float)
        ret[0:3] = self.__RotationMatrixCRF2TRF(epoch).dot(posXyz)
        ret[3:6] = self.__RotationMatrixCRF2TRF(epoch).dot(velXyz)+self.__RotationMatrixDotCRF2TRF(epoch).dot(posXyz)

        return ret

