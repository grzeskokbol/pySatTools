
import datetime
import numpy as np
import libpySat as pySat


class TransformERA:

    def __init__(self,fut1utc):
        self.fut1utc=fut1utc
        self.epochSave = datetime.datetime.now()
        self.rotSave = np.matrix(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
        self.era=0.0
        self.eraDot = 0.0

    def getMatrix_Era(self,epoch: datetime.datetime):
        """
        Functions of Earth Rotation Angle, theta
        Theta is angle bewtween TIO and CIO (along CIP)
        USE UT1 here.
        :param epoch:
        :return:
        """
        if (epoch != self.epochSave):
            mjd = pySat.pySatTime.UT12MJD(epoch)
            ut1utc = self.fut1utc(mjd)
            ut1 = epoch
            ut1+= datetime.timedelta(seconds=float(ut1utc))

            theta, thetadot = self.__getEraDot(epoch, 0.0)

            self.rotSave=pySat.pySatMath.RotationMatrix3DZ(-theta)
            self.epochSave = epoch
            return self.rotSave
        else:
            return self.rotSave

    def __getEraDot(self, epoch: datetime.datetime, ut1utcd):
        """

        :param epoch: epoch in UTC
        :param ut1utcd: temporal differentiation UT1-UTC in seconds/seconds
        :return: theta and thetadot [radians]
        """

        # TODO: Implementation of tidal and libration terms for UT1...
        mjd = pySat.pySatTime.UT12MJD(epoch)
        ut1utc = self.fut1utc(mjd)
        ut1 = epoch + datetime.timedelta(seconds=float(ut1utc))
        # print(ut1)
        # Tu in days
        tu = np.floor(pySat.UT12MJD(ut1)) + pySat.Sod(ut1) / 86400.0 - 51544.5
        # print(tu,pySat.Sod(ut1))
        theta = np.fmod(tu - np.floor(tu) + 0.7790572732640 + 0.00273781191135448 * tu, 1.0) * 2.0 * np.pi
        thetadot = 1.00273781191135448 * 2.0 * np.pi / 86400.0 * (1.0 + ut1utcd)  # rad / utc-sec
        self.era=theta
        self.eraDot=thetadot
        return self.era,self.eraDot

    def getMatrix_EraDot(self, epoch: datetime.datetime):
        """

        :param epoch:
        :return:
        """
        mjd = pySat.pySatTime.UT12MJD(epoch)
        theta,thetadot=self.__getEraDot(epoch,self.fut1utc(mjd))

        return pySat.pySatMath.RotationMatrix3DZ(-theta)*-thetadot

