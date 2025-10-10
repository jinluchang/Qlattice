import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad_vec
from scipy.optimize import curve_fit

class Fit():
    def __init__(self, dt):
        self.dt = dt
    
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
    
    def get_correction(self, t, *argv):
        return self.R0(0,*argv[1:])/self.Rt(np.array([t]), *argv)[0]
    
    def fit_correction(self, t, ts, data, errs, filter_x=None):
        opt = self.get_fit_params(ts,data,errs,filter_x=filter_x)[0]
        return self.get_correction(t, *opt)
    
    def choose_start(self, ts, data, errs, thresh=0.15):
        ts, data, errs = zip(*sorted(zip(ts, data, errs)))
        res_dof = []
        for i in range(len(ts)):
            try:
                info = self.get_fit_params(ts[i:], data[i:], errs[i:], full_output=True)[2]
                res_dof.append(np.linalg.norm(info["fvec"])/(len(ts)-i))
                if(res_dof[-1]<thresh):
                    return ts[i]
            except ValueError:
                print(f"Could not find res/dof below threshold {thresh}. Values were {res_dof} at times {ts}")
                return "err"
    
    def plot_results(self, ts, data, opt, cov, time, dt=1, **plt_args):
        if(data!=[]):
            plt.scatter(ts/dt, data)
        plt.plot(np.divide(ts,dt), self.fit(np.array(ts), *opt), **plt_args)
        print(f"Correction factor: {self.get_correction(time,*opt)}")
    
    def filter_data(self, ts, data, errs, filter_x):
        ts_ = []
        data_ = []
        errs_ = []
        for t in range(len(ts)):
            if filter_x(ts[t]): continue
            ts_.append(ts[t])
            data_.append(data[t])
            errs_.append(errs[t])
        return ts_, data_, errs_

class GaussianFit(Fit):
    def __init__(self, dt, start_time=0):
        self.start_time = start_time
        super().__init__(dt)
    
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
    
    def get_fit_params(self, ts, data, errs, full_output=False, filter_x=None):
        if(filter_x!=None): ts, data, errs = self.filter_data(ts,data,errs,filter_x)
        return curve_fit(self.fit, ts, data, sigma=errs, p0=[1.0, 1.0, 1.0], jac=self.dfit, bounds=((-np.inf,-np.inf,0),(np.inf,np.inf,np.inf)),full_output=full_output)
        #print(np.sqrt(np.diag(cov)))
    
    def get_correction(self, t, *argv):
        argv = (self.start_time,) + argv
        return super().get_correction(t, *argv)
    
    def plot_results(self, ts, data, opt, cov, time, dt=1, **plt_args):
        super().plot_results(ts, data, opt, cov, time, dt, **plt_args)
        try:
            print(f"E_FV: {opt[0]*dt}, E0: {opt[1]*dt}, sigma: {opt[2]*dt}")
        except IndexError:
            print(f"E0: {opt[0]*dt}, sigma: {opt[1]*dt}")
        print(f"Covariance: {np.sqrt(np.diag(cov))}")
    
class GaussianFitNoBounds(GaussianFit):
    def dRt(self, t, start_time, E_FV, *argv):
        # [time, parameter]
        rtn = quad_vec(lambda E: self.dintegrand(E, t, start_time, E_FV, *argv), -E_FV, np.inf)[0]
        return rtn
    
    def fit(self, t, E0, sigma):
        return super().fit(t, 30, E0, sigma)
    
    def dfit(self, t, E0, sigma):
        return super().dfit(t, 30, E0, sigma)[:,1:]
    
    def get_fit_params(self, ts, data, errs, full_output=False, filter_x=None):
        if(filter_x!=None): ts, data, errs = self.filter_data(ts,data,errs,filter_x)
        return curve_fit(self.fit, ts, data, sigma=errs, p0=[1.0, 1.0], jac=self.dfit, bounds=((-np.inf,0),(np.inf,np.inf)), full_output=full_output)
    
    def get_correction(self, t, *argv):
        argv = (30,) + argv
        return super().get_correction(t, *argv)
    
class PowerFit(Fit): 
    def R0(self, E, E_FV, n):
        return np.abs(E+E_FV)**n
    
    def dR0(self, E, E_FV, n):
        return np.column_stack([0, #d_dstart_time
                                n*np.abs(E+E_FV)**(n-1), #d_dE_FV
                                np.log(np.abs(E+E_FV))*np.abs(E+E_FV)**n #d_dn
                               ])
    
    def get_fit_params(self, ts, data, errs, full_output=False, filter_x=None):
        if(filter_x!=None): ts, data, errs = self.filter_data(ts,data,errs,filter_x)
        return curve_fit(self.fit, ts, data, sigma=errs, p0=[1.0, 1.0, 1.0], jac=self.dfit, bounds=((-np.inf,-np.inf,0),(np.inf,np.inf,np.inf)), full_output=full_output)
    
    def plot_results(self, ts, data, opt, cov, time, dt=1.0, **plt_args):
        super().plot_results(ts, data, opt, cov, time, dt, **plt_args)
        print(f"start_time: {opt[0]/dt}, E_FV: {opt[1]*dt}, n: {opt[2]}")
        print(f"Covariance: {np.sqrt(np.diag(cov))}")

class PowerFit2(Fit): 
    def R0(self, E, E_FV, n, a):
        return np.abs(E+E_FV+a)**n
    
    def dR0(self, E, E_FV, n, a):
        return np.column_stack([0, #d_dstart_time
                                n*np.abs(E+E_FV+a)**(n-1), #d_dE_FV
                                np.log(np.abs(E+E_FV+a))*np.abs(E+E_FV+a)**n, #d_dn
                                n*np.abs(E+E_FV+a)**(n-1) #d_da
                               ])
    
    def get_fit_params(self, ts, data, errs, full_output=False, filter_x=None):
        if(filter_x!=None): ts, data, errs = self.filter_data(ts,data,errs,filter_x)
        return curve_fit(self.fit, ts, data, sigma=errs, p0=[1.0, 1.0, 1.0, 0.0], jac=self.dfit, bounds=((-np.inf,-np.inf,0,0),(np.inf,np.inf,np.inf,np.inf)), full_output=full_output)
    
    def plot_results(self, ts, data, opt, cov, time, dt=1.0, **plt_args):
        super().plot_results(ts, data, opt, cov, time, dt, **plt_args)
        print(f"start_time: {opt[0]/dt}, E_FV: {opt[1]*dt}, n: {opt[2]}")
        print(f"Covariance: {np.sqrt(np.diag(cov))}")