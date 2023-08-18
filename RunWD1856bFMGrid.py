import picaso.justdoit as jdi
import picaso.justplotit as jpi
import numpy as np
import astropy.units as u
import pandas as pandas
from myastrotools.reflectx import *

def Run1Model(p, numtangle = 6, numgangle = 6):

    ## Planet:
    #planettype = p['planet_type']
    Tint = p['tint'].item() # Internal Temperature of your Planet in K
    Teq = 163
    radiuse = 10.4 #Rearth
    massj = p['pl_mass'].item() #Mjup
    phase = 0

    r_star = 0.0131 # solar radius
    a_over_rstar = 382.5
    semi_major = (a_over_rstar * r_star)*u.Rsun.to(u.au)

    ## Star:
    star_filename = 'bestfit-JWST-flam-Ang-um-forpicaso.txt'
    
    ## Climate:
    nlevel = int(p['nlevel'].item()) # number of plane-parallel levels in your code
    nofczns = int(p['nofczns'].item()) # number of convective zones initially. Let's not play with this for now.
    nstr_upper = int(p['nstr_upper'].item()) # top most level of guessed convective zone
    nstr_deep = nlevel -2 # this is always the case. Dont change this
    nstr = np.array([0,nstr_upper,nstr_deep,0,0,0]) # initial guess of convective zones
    rfacv = p['rfacv']

    ## Opacities:
    #
    planet_mh = p['mh'].item()
    planet_mh_CtoO = p['cto'].item()

    directory = f'Tint{Tint}-mass{massj}-mh{int(planet_mh)}-co{planet_mh_CtoO}'
    #savefiledirectory = p['output_dir']+directory
    output_dir = '/Volumes/Oy/'
    savefiledirectory = output_dir+directory

    local_ck_path = f'/Volumes/Oy/picaso/reference/kcoeff_2020/'
    #local_ck_path = p['local_ck_path']+'/'

    planet_properties = {
        'tint':Tint, 'Teq':Teq, 'radius':radiuse, 'radius_unit':u.Rearth,
         'mass':massj, 'mass_unit': u.Mjup,
         'gravity': None, 'gravity_unit':None,
        'semi_major':semi_major, 'semi_major_unit': u.AU,
        'mh': planet_mh, 'CtoO':planet_mh_CtoO, 'phase':phase,
        'noTiOVO':p['noTiOVO'].item(), 'planet_mh_str':p['mh_str'],
        'ctostr':p['ctostr'],
        'local_ck_path':local_ck_path, 'num_tangle':numtangle, 'num_gangle':numgangle
    }

    star_properties = {
        'star_filename':star_filename, 'star_radius':r_star
    }

    climate_run_setup = {'climate_pbottom':int(p['p_bottom'].item()),
            'climate_ptop':int(p['p_top'].item()),
            'nlevel':nlevel, 'nofczns':nofczns, 'nstr_upper':nstr_upper,
            'nstr_deep':nstr_deep, 'rfacv':rfacv
    }
    #opa_file = p['opa_file']
   # #opa_file = None
    # wave_range = [float(p['wave_range'].split(',')[0].replace('[','')),
    #           float(p['wave_range'].split(',')[1].replace(' ','').replace(']',''))]
    #wave_range = [0.3,15]
    #spectrum_setup = {'opacity_db':opa_file,
    #                  'wave_range':wave_range,
    #                  'calculation':'reflected', 'R':150
    #                 }

    if p['guess'] == 'guillot':
        use_guillotpt = True

    import os
    os.system('mkdir '+savefiledirectory)
    import time
    k = open(savefiledirectory+'/RunReport.txt','a')
    t = time.localtime()
    outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
    k.write('Started climate run' +outtime + '\n')
    k.close()

    record_terminal_output = True

    try:
        cj = MakeModelCloudFreePlanet(planet_properties, 
                                star_properties,
                                cdict = climate_run_setup,
                                use_guillotpt = use_guillotpt,
                                compute_spectrum = False,
                                specdict = None,
                                savefiledirectory = savefiledirectory,
                                record_terminal_output = record_terminal_output
                    )
        k = open(savefiledirectory+'/RunReport.txt','a')
        t = time.localtime()
        outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
        k.write('Finished ' +outtime)
        k.close()
        
        if record_terminal_output:
            with open(savefiledirectory+'/terminal_output.txt','r') as f:
                z = f.read()
                k = open('WD1856b-RunReport.txt','a')
                t = time.localtime()
                outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
                if 'YAY ! ENDING WITH CONVERGENCE' in z:
                    k.write(savefiledirectory + ' ' +outtime + '  converged  no error \n')
                else:
                    k.write(savefiledirectory + ' ' +outtime + '  FAILED  failed to converge \n')
                k.close()

    except Exception as e:
        k = open('WD1856b-RunReport.txt','a')
        k.write(savefiledirectory + ' ' +outtime + '  FAILED ' +str(e)+' \n')
        k.close()

        k = open(savefiledirectory+'/RunReport.txt','a')
        t = time.localtime()
        outtime = str(t.tm_year)+'-'+str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_hour)+':'+str(t.tm_min)+':'+str(t.tm_sec)
        k.write('Finished ' + outtime + 'with error '+str(e))
        k.close()
        
    return 


def GetP(sheet_id='1BQH36n5O2Kq8iB1ZmM_WNu1RsRVLRpG64ACdKaDYueg', 
             sheet_name='1271251364'):
    url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&gid={sheet_name}"
    p = pd.read_csv(url,dtype={ 'mh_str':np.str_ ,'ctostr':np.str_})
    p = p.dropna(axis=1, how='all')
    for i in range(len(p)):
        try:
            if np.isnan(p[p.columns[0]][0]):
                p = p.drop(i, axis=0)
        except TypeError:
            pass
    return p


def RunGrid(sheet_id='1BQH36n5O2Kq8iB1ZmM_WNu1RsRVLRpG64ACdKaDYueg', 
             sheet_name='1271251364', n_jobs = 6):
    k = open('WD1856b-RunReport.txt','w')
    k.close()
    p = GetP(sheet_id=sheet_id, 
             sheet_name=sheet_name)
    inds = np.array([ 28,  34,  99, 108, 144, 145, 150, 179, 180, 183, 188, 192, 200,
       201, 205, 206, 210, 211, 215, 219, 220, 221, 222, 224, 225, 227,
       230, 231, 232, 235, 236, 239, 240, 242, 243, 251, 260, 261, 262,
       265, 266, 271, 275, 280, 281, 282, 285, 286, 287, 290, 291, 292,
       295, 296, 297, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309,
       310, 311, 312, 313, 315, 316, 320, 321, 325, 326, 335, 340, 341,
       342, 345, 346, 347, 350, 351, 352, 354, 355, 356, 360, 361, 362,
       363, 364, 365, 366, 367, 368, 370, 371, 373, 375, 376, 377, 378,
       380, 381, 382, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393,
       394, 395, 396, 397, 398, 399])
    p = p.loc[inds]
    p = p.reset_index(drop = True)
    import time
    start = time.time()
    
    import picaso.justdoit as jdi
    jdi.Parallel(n_jobs=n_jobs)(jdi.delayed(Run1Model)(p.loc[i]) for i in range(len(p)))
    #Run1Model(p.loc[0])
    stop = time.time()
    import astropy.units as u
    print(stop-start, (stop-start)*u.s.to(u.min))
        
    #return p



RunGrid()
