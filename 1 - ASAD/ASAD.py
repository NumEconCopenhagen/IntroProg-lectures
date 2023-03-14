
from types import SimpleNamespace

import numpy as np

class ASADClass:

    def __init__(self):
        """ initialize the model """

        par = self.par = SimpleNamespace() # parameters
        sim = self.sim = SimpleNamespace() # simulation variables
        datamoms = self.datamoms = SimpleNamespace() # moments in the data
        moms = self.moms = SimpleNamespace() # moments in the model

        # a. externally given parameters
        par.alpha = 0.700 # slope of AD
        par.gamma = 0.075 # slope of SRAS
        par.phi = 0.99 # stickiness in expectations

        # b. parameters to be chosen (here guesses)
        par.delta = 0.80 # AR(1) of demand shock
        par.omega = 0.15 # AR(1) of supply shock
        par.sigma_x = 1.0 # std. of demand shock
        par.sigma_c = 0.2 # st.d of supply shock

        # c. misc paramters
        par.simT = 10_000 # length of simulation

        # d. calculate compound paramters
        self.calc_compound_par()

        # e. simulation
        sim.y_hat = np.zeros(par.simT)
        sim.pi_hat = np.zeros(par.simT)
        sim.z = np.zeros(par.simT)
        sim.x = np.zeros(par.simT)
        sim.s = np.zeros(par.simT)
        sim.c = np.zeros(par.simT)

        # f. data (numbers given in notebook)
        datamoms.std_y = 1.64
        datamoms.std_pi = 0.21
        datamoms.corr_y_pi = 0.31
        datamoms.autocorr_y = 0.84
        datamoms.autocorr_pi = 0.48

    def calc_compound_par(self):
        """ calculates compound parameters """

        par = self.par

        par.a = (1+par.alpha*par.gamma*par.phi)/(1+par.alpha*par.gamma)
        par.beta = 1/(1+par.alpha*par.gamma)

    def simulate(self):
        """ simulate the full model """

        np.random.seed(1917)

        par = self.par
        sim = self.sim

        # a. draw random  shock innovations
        sim.x = np.random.normal(loc=0.0,scale=par.sigma_x,size=par.simT)
        sim.c = np.random.normal(loc=0.0,scale=par.sigma_c,size=par.simT)

        # b. period-by-period
        for t in range(par.simT):

            # i. lagged
            if t == 0:
                z_lag = 0.0
                s_lag = 0.0
                y_hat_lag = 0.0
                pi_hat_lag = 0.0
            else:
                z_lag = sim.z[t-1]
                s_lag = sim.s[t-1]
                y_hat_lag = sim.y_hat[t-1]
                pi_hat_lag = sim.pi_hat[t-1]

            # ii. AR(1) shocks
            z = sim.z[t] = par.delta*z_lag + sim.x[t]
            s = sim.s[t] = par.omega*s_lag + sim.c[t]

            # iii. output and inflation
            sim.y_hat[t] = par.a*y_hat_lag + par.beta*(z-z_lag) \
                - par.alpha*par.beta*s + par.alpha*par.beta*par.phi*s_lag
            sim.pi_hat[t] = par.a*pi_hat_lag + par.gamma*par.beta*z \
                - par.gamma*par.beta*par.phi*z_lag + par.beta*s - par.beta*par.phi*s_lag
            
    def calc_moms(self):
        """ calculate moments """

        # note: same moments as in the data

        sim = self.sim
        moms = self.moms

        moms.std_y = np.std(sim.y_hat)
        moms.std_pi = np.std(sim.pi_hat)
        moms.corr_y_pi = np.corrcoef(sim.y_hat,sim.pi_hat)[0,1]
        moms.autocorr_y = np.corrcoef(sim.y_hat[1:],sim.y_hat[:-1])[0,1]
        moms.autocorr_pi = np.corrcoef(sim.pi_hat[1:],sim.pi_hat[:-1])[0,1]    

    def calc_diff_to_data(self,do_print=False):
        """ calculate difference to data """

        moms = self.moms
        datamoms = self.datamoms

        error = 0.0 # sum of squared differences
        for k in self.datamoms.__dict__.keys():

            diff = datamoms.__dict__[k]-moms.__dict__[k]
            error += diff**2

            if do_print: print(f'{k:12s}| data = {datamoms.__dict__[k]:.4f}, model = {moms.__dict__[k]:.4f}')

        if do_print: print(f'{error = :12.8f}')

        return error