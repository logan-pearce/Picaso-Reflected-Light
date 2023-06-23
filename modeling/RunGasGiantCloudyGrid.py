import numpy as np
#from tools.tools import *
from myastrotools.reflectx import *
import astropy.units as u
import pandas as pandas

def GetP(sheet_id='11u2eirdZcWZdjnKFn3vzqbCtKCodstP-KnoGXC5FdR8', 
             sheet_name='GasGiantsBaseModels'):
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&gid={sheet_name}"
    p = pd.read_csv(url)
    p = p.dropna(axis=1, how='all')
    for i in range(len(p)):
        try:
            if np.isnan(p[p.columns[0]][0]):
                p = p.drop(i, axis=0)
        except TypeError:
            pass
    return p

def Run1CloudyModel(OneCloudp, OneBaseModelp):
    ### Load basemodel:
    p = OneBaseModelp
    grid = p['name']
    # Planet:
    planettype = p['planet_type']
    Teq = ComputeTeq(p['st_teff'], p['rstar']*u.Rsun, p['au']*u.au, 
                     Ab = 0.3, fprime = 1/4) # planet equilibrium temperature 
    radius = p['pl_rad'] #Rjup
    massj = p['pl_mass']
    semi_major = p['au']
    phase = p['phase']
    planet_mh = p['mh']
    planet_mh_CtoO = p['cto']

    # Star:
    T_star = p['st_teff'] # K, star effective temperature
    logg = p['logg'] #logg , cgs
    metal = p['feh'] # metallicity of star
    r_star = p['rstar'] # solar radius

    Teq_str = np.round(Teq, decimals=0)
    directory = f'{grid}-{planettype}-Tstar{int(T_star)}-Rstar{r_star}-Teq{int(Teq_str)}-sep{semi_major}-rad{radius}-mass{massj}-mh{int(planet_mh)}-co{planet_mh_CtoO}-phase{int(phase)}'
    #savefiledirectory = p['output_dir']+directory
    output_dir = ''
    savefiledirectory = output_dir+directory
    
    fsed = OneCloudp['fsed']
    logkzz = np.log10(OneCloudp['kzz'])
    #cloud_filename_prefix = f'cloudy-models/cloudy-fsed{fsed}-kzz{int(logkzz)}'
    cloud_filename_prefix = f'cloudy-fsed{fsed}-kzz{int(logkzz)}'
    #os.system('mkdir '+savefiledirectory+'/cloudy-models/')
    
    
    ### Cloud setup:
    clouds_setup = {'kz':OneCloudp['kzz'], 
                    'fsed':OneCloudp['fsed'], 
                    'mean_mol_weight':OneCloudp['mmw'],
                    'condensates':OneCloudp['condensates'], 
                    'virga_mieff':OneCloudp['virga_mieff']
                   }
    
    ## Spectrum setup
    opa_file = None
    #opa_file = Cloudp['opa_file']
    #R = Cloudp['R']
    R = 150
    wave_range = [float(OneCloudp['wave_range'].split(',')[0].replace('[','')),
              float(OneCloudp['wave_range'].split(',')[1].replace(' ','').replace(']',''))]
    spectrum_setup = {'opacity_db':opa_file,
                      'wave_range':wave_range,
                      'calculation':'reflected', 'R':R
                     }
   
    
    clouds, clouds_added = MakeModelCloudyPlanet(savefiledirectory, 
                                                 clouds_setup, 
                                                 spectrum_setup,
                                                 cloud_filename_prefix,
                                                 calculation = 'planet')


def RunGrid(sheet_id='11u2eirdZcWZdjnKFn3vzqbCtKCodstP-KnoGXC5FdR8', 
             sheet_name='378514385', n_jobs = 3):
    k = open('ReflectXGasGiantRunReport.txt','a')
    k.write('\n ############################### \n')
    k.close()

    CloudyP = GetP(sheet_id=sheet_id, sheet_name='378514385')
    BaseModelP = GetP(sheet_id=sheet_id, sheet_name='0')
    
    import picaso.justdoit as jdi
    #for j in range(len(BaseModelP)):
    for j in range(1):
        jdi.Parallel(n_jobs=n_jobs)(jdi.delayed(Run1CloudyModel)(CloudyP.loc[i], BaseModelP.loc[j]) 
                                    for i in range(3))#range(len(CloudyP)))
        
RunGrid()
        