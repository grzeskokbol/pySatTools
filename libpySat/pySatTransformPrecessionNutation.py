

import datetime
import numpy as np
import libpySat as pySat
from astropy import _erfa as erfa





class TransformPrecessionNutation:

    def __init__(self,fdX,fdY):
        self.fdX=fdX
        self.fdY=fdY
        self.epochSave = datetime.datetime.now()
        self.rotSave = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)


    def __getPoleOffsets(self, epoch: datetime.datetime):
        """
        Gets celestial pole offsets at epoch
        :param epoch:
        :return: dX dY [arcsec]
        """

        mjd = pySat.UTC2MJD(epoch)
        return self.fdX(mjd),self.fdY(mjd)

    def getMatrix_PrecessionNutation(self,epoch: datetime.datetime):
        """
        Returns Precession/Nutation matrix
        :param epoch:
        :return: matrix Q
        """
        if (epoch != self.epochSave):
            Xcio, Ycio, Scio = self.__getCIPandCIOlocator(epoch)
            dX, dY = self.__getPoleOffsets(epoch)
            dX *= np.pi / 180.0 / 3600.0  # from arcsec to radians
            dY *= np.pi / 180.0 / 3600.0
            x = Xcio + dX
            y = Ycio + dY
            s = Scio
            #print(x, y, s, dX, dY, Xcio, Ycio)
            x2 = x * x
            y2 = y * y

            # First matrix of eq(10) / (5.10)
            a = 0.5 + 0.125 * (x2 + y2)
            q1 = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
            q1[0, 0] = 1.0 - a * x2
            q1[0, 1] = -a * x * y
            q1[0, 2] = x
            q1[1, 0] = q1[0, 1]
            q1[1, 1] = 1.0 - a * y2
            q1[1, 2] = y
            q1[2, 0] = -x
            q1[2, 1] = -y
            q1[2, 2] = 1.0 - a * (x2 + y2)
            # print('PrecNutation q1', q1)
            # Second = Rz(s) matrix of eq(10) / (5.10)
            q2 = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
            cc = np.cos(s)
            ss = np.sin(s)
            q2[2, 2] = 1.0
            q2[0, 0] = cc
            q2[1, 1] = cc
            q2[0, 1] = ss
            q2[1, 0] = -ss
            # print('PrecNutation q2', q2)
            self.rotSave=np.matmul(q1, q2)
            self.epochSave = epoch
            return self.rotSave
        else:
            return self.rotSave

    def __getCIPandCIOlocator(self,epoch: datetime.datetime):
        """
        Gets Celestial Intermediate Pole at epoch
        :param epoch:
        :return: X , Y of CIP
        """
        # Gets x,y coords of Celestial Intermediate Pole (CIP) and CIO locator s
        # CIO = Celestial Intermediate Origin
        # Both in GCRS

        mjd = pySat.pySatTime.UTC2MJD(epoch)
        Xcio, Ycio, Scio = erfa.xys06a(2400000.5, mjd)

        return  Xcio, Ycio, Scio
    def getMatrix_PrecessionNutationDot(self,epoch: datetime.datetime):
        """
        Returns Precession/Nutation matrix
        :param epoch:
        :return: matrix Q
        """
        if (epoch != self.epochSave):
            Xcio, Ycio, Scio = self.__getCIPandCIOlocator(epoch)
            dX, dY = self.__getPoleOffsets(epoch)
            dX *= np.pi / 180.0 / 3600.0  # from arcsec to radians
            dY *= np.pi / 180.0 / 3600.0
            x = Xcio + dX
            y = Ycio + dY
            s = Scio
            #print(x, y, s, dX, dY, Xcio, Ycio)
            x2 = x * x
            y2 = y * y

            # First matrix of eq(10) / (5.10)
            a = 0.5 + 0.125 * (x2 + y2)
            q1 = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
            q1[0, 0] = 1.0 - a * x2
            q1[0, 1] = -a * x * y
            q1[0, 2] = x
            q1[1, 0] = q1[0, 1]
            q1[1, 1] = 1.0 - a * y2
            q1[1, 2] = y
            q1[2, 0] = -x
            q1[2, 1] = -y
            q1[2, 2] = 1.0 - a * (x2 + y2)
            # print('PrecNutation q1', q1)
            # Second = Rz(s) matrix of eq(10) / (5.10)
            q2 = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
            cc = np.cos(s)
            ss = np.sin(s)
            q2[2, 2] = 1.0
            q2[0, 0] = cc
            q2[1, 1] = cc
            q2[0, 1] = ss
            q2[1, 0] = -ss
            # print('PrecNutation q2', q2)
            self.rotSave=np.matmul(q1, q2)
            self.epochSave = epoch
            return self.rotSave
        else:
            return self.rotSave