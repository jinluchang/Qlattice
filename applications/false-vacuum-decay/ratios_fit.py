import numpy as np
from scipy.integrate import quad_vec

class Fit():
    def __init__(self, ts, start_time_guess, dt):
        int_range = 1/np.abs(np.min(np.subtract(ts,start_time_guess)))*200
        self.dE = int_range/1000.0
        self.E = np.arange(-int_range,int_range,self.dE)
        self.dt = 0.2
    
    def integrand(self, E, t, start_time, E_FV, *argv):
        #Erange = self.E[self.E>-E_FV]
        # [time, energy]
        #return Erange, np.exp(-(t-start_time)[:,np.newaxis]*Erange[np.newaxis,:])*self.R0(Erange, E_FV, *argv)[np.newaxis,:]
        # [time]
        return np.exp(-(t-start_time)*E)*self.R0(E, E_FV, *argv)
    
    def dintegrand(self, E, t, start_time, E_FV, *argv):
        #Erange = self.E[self.E>-E_FV]
        # [time, energy, parameter]
        #rtn = np.exp(-(t-start_time)[:,np.newaxis,np.newaxis]*Erange[np.newaxis,:,np.newaxis])*self.dR0(Erange, E_FV, *argv)[np.newaxis,:,:]
        #rtn[:,:,0] += Erange[np.newaxis,:]*self.integrand(t, start_time, E_FV, *argv)[1]
        #return Erange, rtn
        # [time, parameter]
        rtn = np.exp(-(t-start_time)*E)[:,np.newaxis]*self.dR0(E, E_FV, *argv)
        rtn[:,0] += E*self.integrand(E, t, start_time, E_FV, *argv)
        return rtn
    
    def Rt(self, t, start_time, E_FV, *argv):
        # [time]
        #return np.sum(self.integrand(t, start_time, E_FV, *argv)[1], axis=1)*self.dE
        return quad_vec(lambda E: self.integrand(E, t, start_time, E_FV, *argv), -E_FV, np.inf)[0]
    
    def dRt(self, t, start_time, E_FV, *argv):
        # [time, parameter]
        #Erange, intgrnd = self.dintegrand(t, start_time, E_FV, *argv)
        #rtn = np.sum(intgrnd, axis=1)*self.dE
        #rtn[:,1] += np.exp(-(t-start_time)*Erange[0])*self.R0(Erange[0], E_FV, *argv)
        #return rtn
        # [time, parameter]
        rtn = quad_vec(lambda E: self.dintegrand(E, t, start_time, E_FV, *argv), -E_FV, np.inf)[0]
        rtn[:,1] += np.exp((t-start_time)*E_FV)*self.R0(-E_FV, E_FV, *argv)
        return rtn

    def fit(self, t, start_time, E_FV, *argv):
        # [time]
        return self.Rt(t-self.dt,start_time,E_FV,*argv)/self.Rt(t,start_time,E_FV,*argv)

    def dfit(self, t, start_time, E_FV, *argv):
        # [time, parameter]
        return (self.dRt(t-self.dt,start_time,E_FV,*argv)*self.Rt(t,start_time,E_FV,*argv)[:,np.newaxis] - self.Rt(t-self.dt,start_time,E_FV,*argv)[:,np.newaxis]*self.dRt(t,start_time,E_FV,*argv)) / self.Rt(t,start_time,E_FV,*argv)[:,np.newaxis]**2

class GaussianFit(Fit):
    def __init__(self, ts, start_time, dt):
        self.start_time = start_time
        super().__init__(ts, start_time, dt)
    
    def R0(self, E, E_FV, E0, sigma):
        return np.exp(-(E-E0)**2/(2*sigma**2))

    def dR0(self, E, E_FV, E0, sigma):
        return np.column_stack([0, # d_dstart_time
                                0, # d_dE_FV
                                2*(E-E0)/(2*sigma**2)*np.exp(-(E-E0)**2/(2*sigma**2)), #d_dE0
                                (E-E0)**2/sigma**3*np.exp(-(E-E0)**2/(2*sigma**2))  #d_dsigma
                               ])
    
    def fit(self, t, E_FV, E0, sigma):
        return super().fit(t, self.start_time, E_FV, E0, sigma)
    
    def dfit(self, t, E_FV, E0, sigma):
        return super().dfit(t, self.start_time, E_FV, E0, sigma)[:,1:]
    
class GaussianFitNoBounds(GaussianFit):
    def dRt(self, t, start_time, E_FV, *argv):
        # [time, parameter]
        rtn = quad_vec(lambda E: self.dintegrand(E, t, start_time, E_FV, *argv), -E_FV, np.inf)[0]
        return rtn
    
    def fit(self, t, E0, sigma):
        return super().fit(t, 1, E0, sigma)
    
    def dfit(self, t, E0, sigma):
        return super().dfit(t, 1, E0, sigma)[:,1:]
    
class PowerFit(Fit): 
    def R0(self, E, E_FV, n):
        return np.abs(E+E_FV)**n
    
    def dR0(self, E, E_FV, n):
        return np.column_stack([0, #d_dstart_time
                                n*np.abs(E+E_FV)**(n-1), #d_dE_FV
                                np.log(np.abs(E+E_FV))*np.abs(E+E_FV)**n #d_dn
                               ])