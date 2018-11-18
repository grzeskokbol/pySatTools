
import datetime
import numpy as np
import libpySat as pySat
from astropy import _erfa as erfa
from scipy.misc import derivative
from scipy import interpolate

class TransformPolarMotion:

    def __init__(self,fxp,fyp):
        self.fxp=fxp
        self.fyp=fyp
        self.epochSave = datetime.datetime.now()
        self.rotSave = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
        self.sprime=0.0

    def __getPolarMotion(self, epoch: datetime.datetime):
        """

        :param epoch:
        :return: polar motion: x,y [mas]

        """
        mjd=pySat.UTC2MJD(epoch)
        return self.fxp(mjd),self.fyp(mjd)

    def __getPolarMotionDot(self, epoch: datetime.datetime):
        """

        :param epoch:
        :return: polar motion: x,y [mas/s]

        """
        mjd=pySat.UTC2MJD(epoch)
        xpdot=derivative(self.fxp,mjd,dx=1e-3,n=1)
        ypdot = derivative(self.fyp, mjd, dx=1e-3,n=1)
        return xpdot,ypdot

    def getMatrix_PolarMotion(self,epoch:datetime.datetime):
        """
        Get the polar motion matrix. Relates ITRF to TIRS.
        :param epoch:
        :return:
        """
        if (epoch !=self.epochSave):
            xp,yp = self.__getPolarMotion(epoch)
            # TODO: Implementation of tidal and libration terms for polar motion...
            xp*=np.pi/180.0/3600.0
            yp*=np.pi/180.0/3600.0
            sp= self.__getTIO(epoch)
            #print(xp,yp,sp)
            rxy= np.matmul(pySat.RotationMatrix3DY(xp),pySat.RotationMatrix3DX(yp))
            rs=pySat.RotationMatrix3DZ(-sp)
            self.rotSave=np.matmul(rs,rxy)
            self.epochSave = epoch
            return self.rotSave
        else:
            return self.rotSave

    def __getTIO(self, epoch:datetime.datetime ):
        """
        Gets the Terrestrial Intermediate Origin (TIO) locator s'
        Terrestrial Intermediate Ref Sys (TIRS) defined by TIO and CIP.
        TIRS related to to CIRS by Earth Rotation Angle
        :param epoch:
        :return:
        """
        mjd = pySat.pySatTime.UTC2MJD(epoch)
        self.sprime=erfa.sp00(2400000.5,mjd)

        return self.sprime

    def getMatrix_PolarMotionDot(self,epoch:datetime.datetime):
        """
        Get the polar motion matrix. Relates ITRF to TIRS.
        :param epoch:
        :return:
        """

        # TODO: Implementation of tidal and libration terms for polar motion...
        xp, yp = self.__getPolarMotion(epoch)
        xpDot,ypDot = self.__getPolarMotionDot(epoch)
        xp *= np.pi / 180.0 / 3600.0
        yp *= np.pi / 180.0 / 3600.0
        xpDot*=np.pi/180.0/3600.0
        ypDot*=np.pi/180.0/3600.0
        spDot = -47.0 / 1.0e6 / 3600.0 / 180.0 * np.pi / 86400.0 / 36525.0
        sp = self.__getTIO(epoch)
        print('Pmotion dot:',xpDot,ypDot,spDot)

        rxy= np.matmul(pySat.RotationMatrix3DY(xp),pySat.RotationMatrix3DX(yp))
        rxyDot = np.matmul(xpDot* pySat.RotationMatrix3DY(xp), pySat.RotationMatrix3DX(yp)) \
                 +np.matmul( pySat.RotationMatrix3DY(xp),ypDot* pySat.RotationMatrix3DX(yp))
        rs=pySat.RotationMatrix3DZ(-sp)
        rsDot=-spDot*pySat.RotationMatrix3DZ(-sp)

        return np.matmul(rsDot,rxy) + np.matmul(rs,rxyDot)
