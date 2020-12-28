version=2.0

import pandas as pd
import os
from scipy.optimize import curve_fit
import numpy as np
import pickle

#FUNCTIONS

#Object to read and store all testrig log files in a directory
class Logfiles(object):
    def __init__(self, working_dir):
        self.version = version
        self.working_dir = working_dir
        fname_list = sorted(os.listdir(working_dir))
        self.vbon_data = pd.DataFrame()
        self.vbon_fnames = []
        self.vbon_indices =[]
        self.vboff_data = pd.DataFrame()
        self.vboff_indices = []
        self.vboff_fnames = []

        for fname in fname_list:
            path = self.working_dir + fname

            parameters = dict(l.strip().split(" = ")  for l in open(path)  if "=" in l)
            testdatetime = pd.to_datetime(parameters["Test Date"] + " "+ parameters["Test Time"])

            x = pd.read_csv(path, skiprows=28, delim_whitespace=True, skipfooter=2, engine='python')

            x.TIME = testdatetime + x.TIME.astype("float")*pd.Timedelta(seconds=1)
            x = x.set_index("TIME")
            x.index.name = ""
            x = x.assign(V = pd.Series(x['ASI/km/h']/3.6).values)

            if 'vbon' in fname:
                self.vbon_indices.append(len(self.vbon_data.index))
                self.vbon_data = self.vbon_data.append(x)
                self.vbon_fnames.append(fname)
            elif 'vboff' in fname:
                self.vboff_indices.append(len(self.vboff_data.index))
                self.vboff_data = self.vboff_data.append(x)
                self.vboff_fnames.append(fname)
            else: print(fname, ": Filename doesn't contain VB state")
        
        
    def calc_CL(self):
        print('Calculating CL vbon....')
        fit = C_fit(self,'CL',True)
        CL_calc = poly_3D(np.array([self.vbon_data['Keel_Angle/deg'],self.vbon_data['V']]),fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7])
        self.vbon_data = self.vbon_data.assign(CL_calc = pd.Series(CL_calc).values)
        self.CLc_vbon = fit
        
        print('Calculating CL vboff....')
        fit = C_fit(self,'CL',False)
        CL_calc = poly_3D(np.array([self.vboff_data['Keel_Angle/deg'],self.vboff_data['V']]),fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7])
        self.vboff_data = self.vboff_data.assign(CL_calc = pd.Series(CL_calc).values)
        self.CLc_vboff = fit
        
        
    def calc_CD(self):
        print('Calculating CD vbon....')
        fit = C_fit(self,'CD',True)
        CD_calc = poly_3D(np.array([self.vbon_data['Keel_Angle/deg'],self.vbon_data['V']]),fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7])
        self.vbon_data = self.vbon_data.assign(CD_calc = pd.Series(CD_calc).values)
        self.CDc_vbon = fit
        
        print('Calculating CD vboff....')
        fit = C_fit(self,'CD',False)
        CD_calc = poly_3D(np.array([self.vboff_data['Keel_Angle/deg'],self.vboff_data['V']]),fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7])
        self.vboff_data = self.vboff_data.assign(CD_calc = pd.Series(CD_calc).values)
        self.CDc_vboff = fit
        
    def calc_CM(self):
        print('Calculating CM vbon....')
        fit = C_fit(self,'CM',True)
        CM_calc = poly_3D(np.array([self.vbon_data['Keel_Angle/deg'],self.vbon_data['V']]),fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7])
        self.vbon_data = self.vbon_data.assign(CM_calc = pd.Series(CM_calc).values)
        self.CMc_vbon = fit
        
        print('Calculating CM vboff....')
        fit = C_fit(self,'CM',False)
        CM_calc = poly_3D(np.array([self.vboff_data['Keel_Angle/deg'],self.vboff_data['V']]),fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7])
        self.vboff_data = self.vboff_data.assign(CM_calc = pd.Series(CM_calc).values)
        self.CMc_vboff = fit
        
    def save_glider (self, fname, name, span, area):
        glider = Glider(self.CLc_vbon,
                        self.CLc_vboff,
                        self.CDc_vbon,
                        self.CDc_vboff,
                        self.CMc_vbon,
                        self.CMc_vboff,
                        name, span, area)
        
        with open(fname + '.pkl', 'wb') as file_out:
            pickle.dump(glider,file_out,-1)
        

#Object to store and process glider data
class Glider(object):
    def __init__(self, CLc_vbon=[],CLc_vboff=[],CDc_vbon=[],CDc_vboff=[],CMc_vbon=[],CMc_vboff=[],name = 'null', span=0,area=0, rho = 1.225):
        self.version = version
        self.CLc_vbon = CLc_vbon
        self.CLc_vboff = CLc_vboff
        self.CDc_vbon = CDc_vbon
        self.CDc_vboff = CDc_vboff
        self.CMc_vbon = CMc_vbon
        self.CMc_vboff = CMc_vboff
        self.name = name
        self.span = span
        self.area = area
        self.rho = rho
    
    def load_glider(self, fname):               
        with open(fname +'.pkl', 'rb') as input:
            loaded=pickle.load(input)
        self.version = loaded.version
        self.CLc_vbon = loaded.CLc_vbon
        self.CLc_vboff = loaded.CLc_vboff
        self.CDc_vbon = loaded.CDc_vbon
        self.CDc_vboff = loaded.CDc_vboff
        self.CMc_vbon = loaded.CMc_vbon
        self.CMc_vboff = loaded.CMc_vboff
        self.name = loaded.name
        self.span = loaded.span
        self.area = loaded.area
            
    def save_glider (self, fname, name):
        glider = Glider(self.CLc_vbon,
                        self.CLc_vboff,
                        self.CDc_vbon,
                        self.CDc_vboff,
                        self.CMc_vbon,
                        self.CMc_vboff,
                        self.name, self.span, self.area)
        
        with open(fname +'.pkl', 'wb') as file_out:
            pickle.dump(glider,file_out,-1)
            
    def glider_specs (self):
        print('Name: ', self.name)
        print('Span: ', self.span)
        print('Area: ', self.area)
        print('CL constants(VG tight): ', self.CLc_vbon)
        print('CL constants(VG loose): ', self.CLc_vboff)
        print('CD constants(VG tight): ', self.CDc_vbon)
        print('CD constants(VG tight): ', self.CDc_vboff)
        print('CM constants(VG tight): ', self.CMc_vbon)
        print('CM constants(VG tight): ', self.CMc_vboff)
        
    
    def calc_Co(self,AoA,V,VG):
        if VG:
            CL_calc = poly_3D(np.array([AoA,V]),self.CLc_vbon[0],self.CLc_vbon[1],self.CLc_vbon[2],self.CLc_vbon[3],self.CLc_vbon[4],self.CLc_vbon[5],self.CLc_vbon[6],self.CLc_vbon[7])
            CD_calc = poly_3D(np.array([AoA,V]),self.CDc_vbon[0],self.CDc_vbon[1],self.CDc_vbon[2],self.CDc_vbon[3],self.CDc_vbon[4],self.CDc_vbon[5],self.CDc_vbon[6],self.CDc_vbon[7])
            CM_calc = poly_3D(np.array([AoA,V]),self.CMc_vbon[0],self.CMc_vbon[1],self.CMc_vbon[2],self.CMc_vbon[3],self.CMc_vbon[4],self.CMc_vbon[5],self.CMc_vbon[6],self.CMc_vbon[7])
        else:
            CL_calc = poly_3D(np.array([AoA,V]),self.CLc_vboff[0],self.CLc_vboff[1],self.CLc_vboff[2],self.CLc_vboff[3],self.CLc_vboff[4],self.CLc_vboff[5],self.CLc_vboff[6],self.CLc_vboff[7])
            CD_calc = poly_3D(np.array([AoA,V]),self.CDc_vboff[0],self.CDc_vboff[1],self.CDc_vboff[2],self.CDc_vboff[3],self.CDc_vboff[4],self.CDc_vboff[5],self.CDc_vboff[6],self.CDc_vboff[7])
            CM_calc = poly_3D(np.array([AoA,V]),self.CMc_vboff[0],self.CMc_vboff[1],self.CMc_vboff[2],self.CMc_vboff[3],self.CMc_vboff[4],self.CMc_vboff[5],self.CMc_vboff[6],self.CMc_vboff[7])
        return CL_calc, CD_calc, CM_calc
            
    
    def calc_LDM(self,AoA,V,VG):
        CL,CD,CM = self.calc_Co(AoA,V,VG)
        q = (self.rho*V**2)/2
        c = self.area/self.span
        L = q * self.area * CL
        D = q * self.area * CD
        M = q * self.area * c * CM
        return L, D, M
            
            
#To fit CL parameters to a logfile
def C_fit(logfiles, C, vb = True):
    if vb:
        data = logfiles.vbon_data
    else:
        data = logfiles.vboff_data
        
    popt, pcov = curve_fit(poly_3D,
                       np.array([data['Keel_Angle/deg'],data['V']]),
                       data[C])                
    sigma = np.sqrt(np.diag(pcov))
    vari = sigma/np.array(popt)
    print('C values: ',popt)
    print('C sigma: ', sigma)
    print('C sigma/values: ', vari)
    print('Mean variance: ',np.mean(vari))
    return popt

#Function to return C as a function of AoA and airspeed
def poly_3D(AoA_V, a,b,c,d, m,n,o,p):
    C = poly_3(AoA_V[0],a,b,c,d)+poly_3(AoA_V[1],m,n,o,p)
    return C

#Polynomial functions to use with curve fitting
def poly_1 (x, lin, const):
    y = lin*x +const
    return y

def poly_2 (x, square, lin, const):
    y = square*x**2 + lin*x +const
    return y

def poly_3 (x, cube, square, lin, const):
    y = cube*x**3 + square*x**2 + lin*x +const
    return y
