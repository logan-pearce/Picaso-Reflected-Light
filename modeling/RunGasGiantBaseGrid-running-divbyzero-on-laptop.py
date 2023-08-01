import numpy as np
#from tools.tools import *
from myastrotools.reflectx import *
import astropy.units as u
import pandas as pandas

def Run1Model(p, num_tangle = 6, num_gangle = 6):

    grid = p['name']
    ## Planet:
    planettype = p['planet_type']
    Tint = p['tint'] # Internal Temperature of your Planet in K
    Teq = ComputeTeq(p['st_teff'], p['rstar']*u.Rsun, p['au']*u.au, 
                     Ab = 0.3, fprime = 1/4) # planet equilibrium temperature 
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
    nlevel = int(p['nlevel']) # number of plane-parallel levels in your code
    nofczns = int(p['nofczns']) # number of convective zones initially. Let's not play with this for now.
    nstr_upper = int(p['nstr_upper']) # top most level of guessed convective zone
    nstr_deep = nlevel -2 # this is always the case. Dont change this
    nstr = np.array([0,nstr_upper,nstr_deep,0,0,0]) # initial guess of convective zones
    rfacv = p['rfacv']

    ## Opacities:
    #
    planet_mh = p['mh']
    planet_mh_CtoO = p['cto']

    Teq_str = np.round(Teq, decimals=0)
    directory = f'{grid}-{planettype}-Tstar{int(T_star)}-Rstar{r_star}-Teq{int(Teq_str)}-sep{semi_major}-rad{radius}-mass{massj}-mh{int(planet_mh)}-co{planet_mh_CtoO}-phase{int(phase)}'
    #savefiledirectory = p['output_dir']+directory
    output_dir = '/Volumes/Oy/ReflectX/ReflectXGasGiantModelGrid/'
    savefiledirectory = output_dir+directory

    local_ck_path = f'/Volumes/Oy/picaso/reference/kcoeff_2020/'
    #local_ck_path = p['local_ck_path']

    planet_properties = {
        'tint':Tint, 'Teq':Teq, 'radius':radius, 'radius_unit':u.Rjup,
         'mass':massj, 'mass_unit': u.Mjup,
         'gravity': None, 'gravity_unit':None,
        'semi_major':semi_major, 'semi_major_unit': u.AU,
        'mh': planet_mh, 'CtoO':planet_mh_CtoO, 'phase':phase, 'num_tangle':num_tangle,
        'num_gangle':num_gangle, 'noTiOVO':p['noTiOVO'], 'planet_mh_str':p['mh_str'],
        'local_ck_path':local_ck_path
    }

    star_properties = {
        'Teff':T_star, 'logg':logg, 'mh':metal, 'radius':r_star
    }

    climate_run_setup = {'climate_pbottom':int(p['p_bottom']),
            'climate_ptop':int(p['p_top']),
            'nlevel':nlevel, 'nofczns':nofczns, 'nstr_upper':nstr_upper,
            'nstr_deep':nstr_deep, 'rfacv':rfacv
    }
    #opa_file = p['opa_file']
    opa_file = '/Volumes/Oy/picaso/reference/references/opacities/all_opacities_0.3_2_R60000.db'
    wave_range = [float(p['wave_range'].split(',')[0].replace('[','')),
              float(p['wave_range'].split(',')[1].replace(' ','').replace(']',''))]
    spectrum_setup = {'opacity_db':opa_file,
                      'wave_range':wave_range,
                      'calculation':'reflected', 'R':p['R']
                     }

    if p['guess'] == 'guillot':
        use_guillotpt = True


    cj = MakeModelCloudFreePlanet(planet_properties, 
                            star_properties,
                            use_guillotpt = True,
                            cdict = climate_run_setup,
                            compute_spectrum = True,
                            specdict = spectrum_setup,
                            savefiledirectory = savefiledirectory
                 )
    import time
    with open(savefiledirectory+'/terminal_output.txt','r') as f:
        z = f.read()
        k = open('ReflectXGasGiantRunReport.txt','a')
        t = time.localtime()
        outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
        if 'YAY ! ENDING WITH CONVERGENCE' in z:
            k.write(savefiledirectory + ' ' +outtime + '  converged \n')
        else:
            k.write(savefiledirectory + ' ' +outtime + '  FAILED \n')
        k.close()
        
    return savefiledirectory, cj


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


def RunGrid(sheet_id='11u2eirdZcWZdjnKFn3vzqbCtKCodstP-KnoGXC5FdR8', 
             sheet_name='GasGiantsBaseModels', n_jobs = 6):
    k = open('ReflectXGasGiantRunReport.txt','w')
    k.close()
    p = GetP(sheet_id=sheet_id, 
             sheet_name=sheet_name)
    locs = [  10,   12,   14,   16,   18,   20,   21,   22,   23,   24,   25,
         26,   27,   28,   29,   30,   30,   31,   31,   32,   32,   33,
         33,   34,   34,   35,   35,   36,   36,   37,   37,   38,   39,
         40,   42,   44,   46,   48,   70,   71,   72,   73,   74,   75,
         76,   77,   78,   79,   80,   81,   82,   83,   84,   85,   86,
         87,   88,   89,  111,  113,  115,  117,  119,  120,  121,  122,
        123,  124,  126,  127,  128,  131,  133,  135,  137,  139,  160,
        161,  162,  164,  166,  168,  170,  172,  174,  176,  178,  184,
        186,  188,  190,  192,  198,  210,  211,  212,  213,  214,  215,
        216,  217,  218,  219,  220,  221,  222,  223,  224,  225,  226,
        227,  228,  229,  230,  232,  234,  236,  238,  290,  292,  293,
        294,  295,  296,  297,  298,  299,  340,  341,  342,  343,  344,
        345,  346,  348,  349,  430,  432,  434,  436,  438,  440,  441,
        443,  447,  449,  480,  484,  486,  488,  941,  943,  945,  947,
        990,  991,  992,  994,  995,  996,  997,  998,  999, 1004, 1006,
       1051, 1054, 1055, 1059, 1101, 1105, 1107, 1150, 1152, 1156, 1157,
       1158, 1208, 1209]
    p = p.loc[locs]
    p = p.reset_index(drop=True)
    import picaso.justdoit as jdi
    jdi.Parallel(n_jobs=n_jobs)(jdi.delayed(Run1Model)(p.loc[i]) for i in range(len(p)))
        
    #return p



RunGrid()