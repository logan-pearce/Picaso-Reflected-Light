import numpy as np
#from tools.reflectx import *
from tools.reflectx import *
import astropy.units as u
import pandas as pd

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
    directory = f'{grid}-{planettype}-Tstar{int(T_star)}-Rstar{r_star}-Teq{int(Teq_str)}-sep{semi_major}-rad{radius}-mass{massj}-mh{planet_mh}-co{planet_mh_CtoO}-phase{int(phase)}'
    savefiledirectory = p['output_dir']+directory
    #output_dir = ''
    #savefiledirectory = output_dir+directory

    # Skip ones already completed:
    # try:
    #     with open('ReflectXGasGiantRunReport.txt','a') as f:
    #         z = f.read().splitlines()
    #     for zz in z:
    #         if directory in zz:
    #             return
    # except:
    #     pass

    import os
    # Make directory to store run results:
    #os.system('mkdir '+savefiledirectory)

    import time
    os.system('touch '+savefiledirectory+'/RunReport.txt')
    k = open(savefiledirectory+'/RunReport.txt','w')
    t = time.localtime()
    outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
    k.write('Starting ' +outtime + '\n')
    k.close()


    #local_ck_path = f'/Volumes/Oy/picaso/reference/kcoeff_2020/'
    local_ck_path = p['local_ck_path']

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

    ### V2 change:
    if int(p['p_bottom']) == 2:
        # No need to rerun the ok ones from V1.
        k = open('ReflectXGasGiantRunReport.txt','a')
        k.write(savefiledirectory + '  Use V1 PT Profile bc it was fine \n')
        k.close()
        return
    if int(p['p_bottom']) == 3:
        # set higher pbottom to 300 instead of 1000 cause that's too much:
        p_bottom = np.log10(300)
    else:
        p_bottom = float(p['p_bottom'])

    climate_run_setup = {'climate_pbottom':p_bottom,
            'climate_ptop':int(p['p_top']),
            'nlevel':nlevel, 'nofczns':nofczns, 'nstr_upper':nstr_upper,
            'nstr_deep':nstr_deep, 'rfacv':rfacv
    }
    opa_file = p['opa_file']
    #opa_file = None
    R = p['R']
    wave_range = [float(p['wave_range'].split(',')[0].replace('[','')),
            float(p['wave_range'].split(',')[1].replace(' ','').replace(']',''))]
    spectrum_setup = {'opacity_db':opa_file,
                    'wave_range':wave_range,
                    'calculation':'reflected', 'R':R
                    }

    if p['guess'] == 'guillot':
        use_guillotpt = True
    
    import time
    k = open(savefiledirectory+'/RunReport.txt','a')
    t = time.localtime()
    outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
    k.write('Started climate run' +outtime + '\n')
    k.close()

    try:
        pl, noclouds, recommended = MakeModelCloudFreePlanet(planet_properties, 
                                star_properties,
                                cdict = climate_run_setup,
                                use_guillotpt = use_guillotpt,
                                #compute_spectrum = True,
                                # V2 change:
                                compute_spectrum = False,
                                specdict = spectrum_setup,
                                savefiledirectory = savefiledirectory,
                                record_terminal_output = True
                    )
        k = open(savefiledirectory+'/RunReport.txt','a')
        t = time.localtime()
        outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
        k.write('Finished ' +outtime)
        k.close()
        
        with open(savefiledirectory+'/terminal_output.txt','r') as f:
            z = f.read()
            k = open('ReflectXGasGiantRunReport.txt','a')
            t = time.localtime()
            outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
            if 'YAY ! ENDING WITH CONVERGENCE' in z:
                k.write(savefiledirectory + ' ' +outtime + '  converged  no error \n')
            else:
                k.write(savefiledirectory + ' ' +outtime + '  FAILED  failed to converge \n')
            k.close()

        GenerateInitialReadMe(savefiledirectory,planettype,T_star,r_star,Teq_str,semi_major,radius,massj,
                   planet_mh,planet_mh_CtoO,phase,planet_properties,
                   star_properties,climate_run_setup,spectrum_setup,recommended)

    except Exception as e:
        k = open('ReflectXGasGiantRunReport.txt','a')
        k.write(savefiledirectory + ' ' +outtime + '  FAILED ' +str(e)+' \n')
        k.close()

        k = open(savefiledirectory+'/RunReport.txt','a')
        t = time.localtime()
        outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
        k.write('Finished ' + outtime + 'with error '+str(e))
        k.close()
        
    return 


def GetP(sheet_id='1rD4aVpD57SQuPR8f2cqb2IJT67KlOg9NNnOu0lbUqGg', 
             sheet_name=0):
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


def RunGrid(sheet_id='1rD4aVpD57SQuPR8f2cqb2IJT67KlOg9NNnOu0lbUqGg', 
             sheet_name=0, n_jobs = 10):
    k = open('ReflectXGasGiantRunReport.txt','a')
    k.close()
    p = GetP(sheet_id=sheet_id, 
             sheet_name=sheet_name)
    #p = p.loc[:0]
    # import picaso.justdoit as jdi
    # jdi.Parallel(n_jobs=n_jobs)(jdi.delayed(Run1Model)(p.loc[i]) for i in range(len(p)))

    from multiprocessing import Pool
    n_cpu = n_jobs
    listofthingstodo = [p.loc[i] for i in range(len(p))]
    with Pool(n_cpu) as pool:
        result = pool.map(Run1Model, listofthingstodo)
    #Run1Model(p)
        
    #return p



RunGrid()