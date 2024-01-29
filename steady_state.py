# find steady state

import time
import numpy as np

from EconModel import jit

from consav.grids import equilogspace
from consav.markov import tauchen, find_ergodic
from consav.misc import elapsed
from scipy import optimize

import household_problem

def prepare_hh_ss(model):
    """ prepare the household block for finding the steady state """
    
    par = model.par
    ss = model.ss
    
    ############
    # 1. grids #
    ############
    
    par.a_grid[:] = equilogspace(0.0,par.a_max,par.Na)
    
    par.beta_shares[:] = np.array([par.HtM_share,1-par.HtM_share-par.PIH_share,par.PIH_share])
    
    ###########################
    # 2. initial distribution #
    ###########################
    
    for i_fix in range(par.Nfix):
        ss.Dbeg[i_fix,:,:] = 0.0      
        ss.Dbeg[i_fix,0,0] = par.beta_shares[i_fix]
    
    ################################################
    # 3. initial guess for intertemporal variables #
    ################################################
    
    model.set_hh_initial_guess()
    
def find_ss(model,do_print=False):
    
    par = model.par
    ss = model.ss 

    ss.w = par.w_ss
    
    ss.delta = delta_ss = 0.02748598211864653
    ss.lambda_u_s = 0.3056686310719316
    
    ss.UI_ratio_high = par.UI_ratio_high
    ss.UI_duration = par.UI_duration
    ss.tau = par.tau       
    
    ###### HOUSEHOLD PROBLEM ######
                
    model.solve_hh_ss()
    model.simulate_hh_ss()
    
    ###### POST HOUSEHOLD PROBLEM ######
    
    # Step 1: 
    ss.u = ss.U_POL_hh

    ss.S = ss.S_hh = np.sum(ss.Dbeg*ss.s)
    ss.ut = (1-ss.u)*ss.delta + ss.u

  