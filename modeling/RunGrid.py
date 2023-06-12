import numpy as np
from tools.tools import *
import astropy.units as u
import pandas as pandas

def Run1Model(p):

    grid = p['name']
    ## Planet:
    planettype = p['planet_type']
    Tint = p['tint'] # Internal Temperature of your Planet in K
    Teq = ComputeTeq(p['st_teff'], p['rstar'], p['au'], Ab = 0.3, fprime = 1/4) # planet equilibrium temperature 
    radius = p['pl_rad'] #Rjup
    massj = p['pl_mass']
    semi_major = p['au']
    phase = p['phase']

    ## Star:
    T_star = p['st_teff'] # K, star effective temperature
    logg = p['logg'] #logg , cgs
    metal = p['feh'] # metallicity of star
    r_star = p['rstar'] # solar radius


    ## Climate:
    nlevel = p['nlevel'] # number of plane-parallel levels in your code
    nofczns = p['nofczns'] # number of convective zones initially. Let's not play with this for now.
    nstr_upper = p['nstr_upper'] # top most level of guessed convective zone
    nstr_deep = nlevel -2 # this is always the case. Dont change this
    nstr = np.array([0,nstr_upper,nstr_deep,0,0,0]) # initial guess of convective zones
    rfacv = p['rfacv']

    ## Opacities:
    #
    planet_mh = p['mh']
    planet_mh_str = p['mh_str']
    #
    planet_mh_CtoO = p['cto']

    directory = f'{grid}-{planettype}-Tstar{T_star}-Rstar{r_star}-Teq{Teq}-sep{semi_major}-rad{radius}-mass{massj}-mh{planet_mh}-co{planet_mh_CtoO}-phase{phase}'
    savefiledirectory = p['output_dir']+directory


    if p['noTiOVO']:
        ck_db_name = p['local_ck_path'] + f'sonora_2020_feh{planet_mh_str}_co_{planet_mh_CtoO}_noTiOVO.data.196'
    else:
        ck_db_name = p['local_ck_path'] + f'sonora_2020_feh{planet_mh_str}_co_{planet_mh_CtoO}.data.196'
        
        
    planet_properties = {
        'tint':Tint, 'Teq':Teq, 'radius':radius, 'radius_unit':u.Rjup,
         'mass':massj, 'mass_unit': u.Mjup,
         'gravity': None, 'gravity_unit':None,
        'semi_major':semi_major, 'semi_major_unit': u.AU,
        'mh': planet_mh_str, 'CtoO':planet_mh_CtoO, 'phase':phase
    }

    star_properties = {
        'Teff':T_star, 'logg':logg, 'mh':metal, 'radius':r_star
    }

    climate_run_setup = {
        'nlevel':nlevel, 'nofczns':nofczns, 'nstr_upper':nstr_upper,
        'nstr_deep':nstr_deep, 'rfacv':rfacv
    }

    if p['guess'] == 'guillot':
        use_guillotpt = True


    cj = MakeModelCloudFreePlanet(planet_properties, 
                            star_properties, ck_db_name,
                            use_guillotpt = True,
                            cdict = climate_run_setup,
                            climate_pbottom = p['p_bottom'],
                            climate_ptop = p['p_top'], 
                            savemodel = True,
                            savefiledirectory = savefiledirectory
                 )

def RunGrid(sheet_id='11u2eirdZcWZdjnKFn3vzqbCtKCodstP-KnoGXC5FdR8', 
             sheet_name='GasGiantsBaseModels'):

    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&gid={sheet_name}"
    p = pd.read_csv(url)
    p = p.dropna(axis=1, how='all')
    print(p.columns)



RunGrid()


