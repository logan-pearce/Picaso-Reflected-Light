import numpy as np
import astropy.units as u
import astropy.constants as c
import os
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

def ComputeTeq(StarTeff, StarRad, sep, Ab = 0.3, fprime = 1/4):
    ''' from Seager 2016 Exoplanet Atmospheres eqn 3.9
    https://books.google.com/books?id=XpaYJD7IE20C
    '''
    StarRad = StarRad.to(u.km)
    sep = sep.to(u.km)
    return (StarTeff * np.sqrt(StarRad/sep) * ((fprime * (1 - Ab))**(1/4))).value

def GetPhaseAngle(MeanAnom, ecc, inc, argp):
    ''' Function for returning observed phase angle given orbital elements
    Args:
        MeanAnom (flt): Mean anomly in radians, where MeanAnom = orbit fraction*2pi, or M=2pi * time/Period
        ecc (flt): eccentricity, defined on [0,1)
        inc (flt): inclination in degrees, where inc = 90 is edge on, inc = 0 or 180 is face on orbit
        argp (flt): argument of periastron in degrees, defined on [0,360)
        
    Returns:
        flt: phase angle in degrees
    Written by Logan Pearce, 2023
    '''
    import numpy as np
    from myastrotools.tools import danby_solve, eccentricity_anomaly
    inc = np.radians(inc)
    argp = np.radians(argp)
    EccAnom = danby_solve(eccentricity_anomaly, MeanAnom, ecc, 0.001, maxnum=50)
    TrueAnom = 2*np.arctan( np.sqrt( (1+ecc)/(1-ecc) ) * np.tan(EccAnom/2) )
    Alpha = np.arccos( np.sin(inc) * np.sin(TrueAnom + argp) )
    
    return np.degrees(Alpha)



def MakeModelPlanet(pdict, sdict, opacity_db,
                calculation = "planet",
                use_guillotpt = True,
                user_supplied_ptprofile = None,
                compute_climate = True,
                cdict = None,
                climate_pbottom = 2,
                climate_ptop = -6, 
                spectrum_wavelength_range = [0.5,1.8],
                spectrum_calculation = 'reflected',
                spectrum_resolution = 150,
                add_clouds = True,
                clouddict = None,
                molecules = None,
                savemodel = False,
                savefilename = None
             ):
    
    ''' Wrapper for PICASO functions for building a planet model
    Args:
        pdict (dict): dictionary of planet parameter inputs
        sdict (dict): dictionary of star parameter inputs
        opacity_db (jdi.opannection object)
        calculation (str): picaso input for object, "planet" or "brown dwarf"
        use_guillotpt (bool): if True, use Guillot PT approximation. Else user must supply initial PT profile
        user_supplied_ptprofile (df): user supplied pt profile for picaso
        comput_climate (bool): if true use picaso to compute plnet climate
        cdict (dict): dictionary of climate run setup params
        climate_pbottom (flt): log(pressure) at bottom of climate calc
        climate_ptop (flt): log(pressure) at top of climate calc
        spectrum_wavelength_range (list): range in um of wavelengths to compute spectrum
        spectrum_calculation (str): type of spectrum to calculate
        spectrum_resolution (flt): what R to compute the spectrum
        add_clouds (bool): If true, add clouds to model
        clouddict (dict): dictionary of cloud parameters
        molecules (list): list of molecules to compute cloud properties. If None, use virga recommended mols
        savemodel (bool): if true, save the model using the xarray method in picaso
        savefilename (str): filename and path for the model to be saved.
    Returns:
        pl: picaso planet model inputs
        noclouds: picaso object after climate run before clouds
        w_noclouds, f_noclouds: wavelength and flux arrays for noclouds spectrum sampled at spectrum_resolution
        clouds_added: virga output from adding clouds
        clouds_spectrum: spectrum after adding clouds computed from spectrum_calculation
        w_clouds, f_clouds: cloudy spectrum sampled at spectrum_resolution

    '''
    import warnings
    warnings.filterwarnings('ignore')
    import picaso.justdoit as jdi

    add_output={
            'author':"Logan Pearce",
            'contact' : "loganpearce1@arizona.edu",
            'code' : "picaso, virga",
            'planet_params':pdict,
            'stellar_params':sdict,
            'orbit_params':{'sma':pdict['semi_major']}
            }
    
    # initialize model:
    pl = jdi.inputs(calculation= calculation, climate = compute_climate) # start a calculation
    # set up planet:
    pl.effective_temp(pdict['tint']) # input effective temperature
    # add gravity:
    if not pdict['gravity']:
        pl.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        pl.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
        
    # set up star:
    pl.star(opacity_db, temp = sdict['Teff'], metal = sdict['mh'], logg = sdict['logg'], 
            radius = sdict['radius'], radius_unit = u.R_sun, 
            semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'], database = 'phoenix')
    
    # climate run
    if use_guillotpt:
        pt = pl.guillot_pt(pdict['Teq'], nlevel=cdict['nlevel'], T_int = pdict['tint'], 
                              p_bottom=climate_pbottom, p_top=climate_ptop)
    else:
        pt = user_supplied_ptprofile

    if compute_climate:
        # initial PT profile guess:
        temp_guess = pt['temperature'].values 
        press_guess = pt['pressure'].values
        # Input climate params:
        nstr = np.array([0,cdict['nstr_upper'],cdict['nstr_deep'],0,0,0]) # initial guess of convective zones
        pl.inputs_climate(temp_guess= temp_guess, pressure= press_guess, 
                      nstr = nstr, nofczns = cdict['nofczns'] , rfacv = cdict['rfacv'])
        print('starting climate run')
        # Compute climate:
        noclouds = pl.climate(opacity_db, save_all_profiles=True, with_spec=True)
    else:
        noclouds = pl.copy()
        
    # make a new object for computing the new spectrum:
    opa_mon = jdi.opannection(wave_range=spectrum_wavelength_range)
    noclouds_spec = jdi.inputs(calculation="planet") # start a calculation
    noclouds_spec.phase_angle(0)
    # add gravity:
    if not pdict['gravity']:
        noclouds_spec.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        noclouds_spec.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
        # add same star:
    noclouds_spec.star(opa_mon, temp = sdict['Teff'], metal = sdict['mh'], logg = sdict['logg'], 
        radius = sdict['radius'], radius_unit=u.R_sun, 
        semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'])
    # add new atmosphere computer by climate run:
    noclouds_spec.atmosphere(df=noclouds['ptchem_df'])
    # compute spectrum:
    noclouds_spec_spectrum = noclouds_spec.spectrum(opa_mon, 
                                                    calculation=spectrum_calculation, 
                                                    full_output=True)
    w_noclouds, f_noclouds = jdi.mean_regrid(noclouds_spec_spectrum['wavenumber'],
                          noclouds_spec_spectrum['fpfs_reflected'], R=spectrum_resolution)
    if not add_clouds:
        # if no clouds, save model and finish computation:
        if savemodel:
            preserve = jdi.output_xarray(
                noclouds_spec_spectrum,
                noclouds_spec,
                add_output = add_output,
                savefile=savefilename
                )
        return pl, noclouds, noclouds_spec_spectrum, w_noclouds, f_noclouds
    else:
        # else add clouds:
        print('adding clouds')
        from virga import justdoit as vj
        # pressure temp profile from climate run:
        temperature = noclouds['temperature']
        pressure = noclouds['pressure']
        #metallicity = pdict['mh'] #atmospheric metallicity relative to Solar
        metallicity_TEMP = 0
        # got molecules for cloud species:
        if not molecules:
            # if no user-supplied molecules:
            #get virga recommendation for which gases to run
            # metallicitiy must be in NOT log units
            recommended = vj.recommend_gas(pressure, temperature, 10**(metallicity_TEMP), clouddict['mean_mol_weight'],
                            #Turn on plotting
                             plot=False)
            mols = recommended
            print('virga recommended gas species:',recommended)
        else:
            # otherwise use user supplied mols:
            mols = molecules
        # add atmosphere from climate run:
        atm = noclouds['ptchem_df']
        # add kzz:
        atm['kz'] = [clouddict['kz']]*atm.shape[0]
        
        # Just set all this up again:
        clouds = jdi.inputs(calculation=calculation) # start a calculation
        clouds.phase_angle(0)
        # add gravity:
        if not pdict['gravity']:
            clouds.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
                mass = pdict['mass'], mass_unit=pdict['mass_unit'])
        else:
            clouds.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
        # add star:
        clouds.star(opa_mon, temp = sdict['Teff'], metal = sdict['mh'], logg = sdict['logg'], 
            radius = sdict['radius'], radius_unit=u.R_sun, 
            semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'], database = 'phoenix')
        # add atmosphere from climate run with kzz:
        clouds.atmosphere(df=atm)
        # get clouds from reference:
        directory ='/Volumes/Oy/virga/virga/reference/RefIndexFiles'
        clouds_added = clouds.virga(mols,directory, fsed=clouddict['fsed'], mh=10**(metallicity_TEMP),
                 mmw = clouddict['mean_mol_weight'], full_output=True)
        # compute spectrum:
        clouds_spectrum = clouds.spectrum(opa_mon, 
                                                calculation=spectrum_calculation, 
                                                full_output=True)
        w_clouds,f_clouds = jdi.mean_regrid(clouds_spectrum['wavenumber'],
                      clouds_spectrum['fpfs_reflected'], R=spectrum_resolution)
        # save:
        if savemodel:
            preserve = jdi.output_xarray(
                clouds_spectrum,
                clouds,
                add_output = add_output,
                savefile=savefilename)
            import pickle
            pickle.dump([pl, noclouds, w_noclouds, f_noclouds, clouds, clouds_added, mols, clouds_spectrum, w_clouds, f_clouds],
                        open(savefilename.replace('.nc','.pkl'),'wb'))
            
        
        return pl, noclouds, w_noclouds, f_noclouds, clouds, clouds_added, mols, clouds_spectrum, w_clouds, f_clouds


##### Broken into modules:
def MakeModelCloudFreePlanet(pdict, sdict,
                calculation = "planet",
                use_guillotpt = True,
                user_supplied_ptprofile = None,
                cdict = None,
                climate_pbottom = 2,
                climate_ptop = -6,
                savemodel = False,
                savefiledirectory = None
             ):
    
    ''' Wrapper for PICASO functions for building a planet model
    Args:
        pdict (dict): dictionary of planet parameter inputs
        sdict (dict): dictionary of star parameter inputs
        opacity_db (jdi.opannection object)
        calculation (str): picaso input for object, "planet" or "brown dwarf"
        use_guillotpt (bool): if True, use Guillot PT approximation. Else user must supply initial PT profile
        user_supplied_ptprofile (df): user supplied pt profile for picaso
        cdict (dict): dictionary of climate run setup params
        climate_pbottom (flt): log(pressure) at bottom of climate calc
        climate_ptop (flt): log(pressure) at top of climate calc
        molecules (list): list of molecules to compute cloud properties. If None, use virga recommended mols
        savemodel (bool): if true, save the model using the xarray method in picaso
        savefilename (str): filename and path for the model to be saved.
    Returns:
        pl: picaso planet model inputs
        noclouds: picaso object after climate run before clouds
    '''
    import warnings
    warnings.filterwarnings('ignore')
    import picaso.justdoit as jdi
    
    import sys
    import os
    os.system('mkdir '+savefiledirectory)
    

    f = open(savefiledirectory+"/terminal_output.txt", 'w')
    sys.stdout = f

    add_output={
            'author':"Logan Pearce",
            'contact' : "loganpearce1@arizona.edu",
            'code' : "picaso, virga",
            'planet_params':pdict,
            'stellar_params':sdict,
            'orbit_params':{'sma':pdict['semi_major']}
            }
    
    # retrieve opacity correlated k-tables database:
    PlanetMH = pdict['mh']
    PlanetCO = pdict['CtoO']
    ck_db = f'/Volumes/Oy/picaso/reference/kcoeff_2020/sonora_2020_feh{PlanetMH}_co_{PlanetCO}.data.196'
    opacity_ck = jdi.opannection(ck_db=ck_db)
    
    # initialize model:
    pl = jdi.inputs(calculation= calculation, climate = True)
    
    # set up planet:
    # input effective temperature
    pl.effective_temp(pdict['tint']) 
    # add gravity:
    if not pdict['gravity']:
        pl.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        pl.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
        
    # set up star:
    pl.star(opacity_ck, temp = sdict['Teff'], metal = sdict['mh'], logg = sdict['logg'], 
            radius = sdict['radius'], radius_unit = u.R_sun, 
            semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'], database = 'phoenix')
    
    # climate run
    if use_guillotpt:
        pt = pl.guillot_pt(pdict['Teq'], nlevel=cdict['nlevel'], T_int = pdict['tint'], 
                              p_bottom=climate_pbottom, p_top=climate_ptop)
    else:
        pt = user_supplied_ptprofile

    # initial PT profile guess:
    temp_guess = pt['temperature'].values 
    press_guess = pt['pressure'].values
    # Input climate params:
    nstr = np.array([0,cdict['nstr_upper'],cdict['nstr_deep'],0,0,0]) # initial guess of convective zones
    pl.inputs_climate(temp_guess= temp_guess, pressure= press_guess, 
                  nstr = nstr, nofczns = cdict['nofczns'] , rfacv = cdict['rfacv'])
    print('starting climate run')
    # Compute climate:
    noclouds = pl.climate(opacity_ck, save_all_profiles=True, with_spec=True)

    from virga import justdoit as vj
    # pressure temp profile from climate run:
    temperature = noclouds['temperature']
    pressure = noclouds['pressure']
    #metallicity = pdict['mh'] #atmospheric metallicity relative to Solar
    metallicity_TEMP = 0
    mmw = 2.2
    # got molecules for cloud species:
    recommended, fig = vj.recommend_gas(pressure, temperature, 10**(metallicity_TEMP), 
                                   mmw, plot=True, outputplot = True)
    
    
    from bokeh.plotting import figure, output_file, save
    output_file(savefiledirectory+"/recomended-gasses.html")
    save(fig)
    
    import pickle
    pickle.dump([pl, noclouds], open(savefiledirectory+'/cloud-free-model.pkl','wb'))
    pickle.dump([pdict, sdict, cdict], open(savefiledirectory+'/cloud-free-model-inputs.pkl','wb'))
    
    f.close()
    
    
    return pl, noclouds
    

def MakeModelCloudyPlanet(savefiledirectory, clouddict,
                          calculation = 'planet', 
                         molecules = None):
    
    ''' Wrapper for PICASO functions for building a planet model
    Args:
        savefiledirectory (str): directory containing picaso cloud-free model base case.
        clouddict (dict): dictionary of cloud parameter inputs
        calculation (str): picaso input for object, "planet" or "brown dwarf"
        molecules (list): list of molecules to compute cloud properties. If None, use virga recommended mols
    Returns:
        clouds: picaso planet model inputs
        clouds_added: virga cloud run output clouds
    '''
    # import climate run output:
    import pickle
    import picaso.justdoit as jdi
    
    import sys
    f = open(savefiledirectory+"/terminal_output.txt", 'a')
    sys.stdout = f
    print()
    
    pl, noclouds = pickle.load(open(savefiledirectory+'/cloud-free-model.pkl','rb'))
    pdict, sdict, cdict = pickle.load(open(savefiledirectory+'/cloud-free-model-inputs.pkl','rb'))
    
    opa_mon = jdi.opannection()

    from virga import justdoit as vj
    # pressure temp profile from climate run:
    temperature = noclouds['temperature']
    pressure = noclouds['pressure']
    #metallicity = pdict['mh'] #atmospheric metallicity relative to Solar
    metallicity_TEMP = 0
    # got molecules for cloud species:
    if not molecules:
        # if no user-supplied molecules:
        #get virga recommendation for which gases to run
        # metallicitiy must be in NOT log units
        recommended = vj.recommend_gas(pressure, temperature, 10**(metallicity_TEMP), 
                                       clouddict['mean_mol_weight'],
                        #Turn on plotting
                         plot=False)
        mols = recommended
        print('using virga recommended gas species:',recommended)
    else:
        # otherwise use user supplied mols:
        mols = molecules
    # add atmosphere from climate run:
    atm = noclouds['ptchem_df']
    # add kzz:
    atm['kz'] = [clouddict['kz']]*atm.shape[0]

    # Just set all this up again:
    clouds = jdi.inputs(calculation=calculation) # start a calculation
    clouds.phase_angle(0)
    # add gravity:
    if not pdict['gravity']:
        clouds.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        clouds.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
    # add star:
    clouds.star(opa_mon, temp = sdict['Teff'], metal = sdict['mh'], logg = sdict['logg'], 
        radius = sdict['radius'], radius_unit=u.R_sun, 
        semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'])
    # add atmosphere from climate run with kzz:
    clouds.atmosphere(df=atm)
    # get clouds from reference:
    directory ='/Volumes/Oy/virga/virga/reference/RefIndexFiles'
    
    clouds_added = clouds.virga(mols,directory, fsed=clouddict['fsed'], mh=clouddict['mh'],
                         mmw = clouddict['mean_mol_weight'], full_output=True)

    clouddict.update({'condensates':mols})

    pickle.dump(clouds,
                    open(savefiledirectory+'/cloudy-model.pkl','wb'))
    pickle.dump([pdict, sdict, cdict, clouddict],open(savefiledirectory+'/cloudy-model-inputs.pkl','wb'))
    
    f.close()

    return clouds, clouds_added


def MakeModelCloudyAndCloudFreeSpectra(savefiledirectory,
                            spectrum_wavelength_range = [0.5,1.8],
                            spectrum_calculation = 'reflected',
                            spectrum_resolution = 150,
                            calculation = "planet",
                            plot_albedo = False
                                      ):
    
    ''' Wrapper for PICASO functions for building a planet model
    Args:
        savefiledirectory (str): directory containing picaso cloud-free model base case.
        spectrum_wavelength_range (list): range in um of wavelengths to compute spectrum
        spectrum_calculation (str): type of spectrum to calculate
        spectrum_resolution (flt): what R to compute the spectrum
        calculation (str): picaso input for object, "planet" or "brown dwarf"
        plot_albedo (bool): if True, return spectrum in albedo, otherwise return planet/star flux ratio
    Returns:
        clouds: picaso planet model inputs
        clouds_added: virga cloud run output clouds
    '''
    import pickle
    import picaso.justdoit as jdi
    import matplotlib.pyplot as plt

    opa_mon = jdi.opannection(wave_range=spectrum_wavelength_range)
    
    ### Cloud-free spectrum:
    pl, noclouds = pickle.load(open(savefiledirectory+'/cloud-free-model.pkl','rb'))
    pdict, sdict, cdict = pickle.load(open(savefiledirectory+'/cloud-free-model-inputs.pkl','rb'))
    noclouds_spec = jdi.inputs(calculation="planet") # start a calculation
    noclouds_spec.phase_angle(0)
    # add gravity:
    if not pdict['gravity']:
        noclouds_spec.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        noclouds_spec.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
        # add same star:
    noclouds_spec.star(opa_mon, temp = sdict['Teff'], metal = sdict['mh'], logg = sdict['logg'], 
        radius = sdict['radius'], radius_unit=u.R_sun, 
        semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'])
    # add new atmosphere computer by climate run:
    noclouds_spec.atmosphere(df=noclouds['ptchem_df'])
    # compute spectrum:
    noclouds_spec_spectrum = noclouds_spec.spectrum(opa_mon, 
                                                    calculation=spectrum_calculation, 
                                                    full_output=True)
    if plot_albedo:
        w_noclouds, f_noclouds = jdi.mean_regrid(noclouds_spec_spectrum['wavenumber'],
                              noclouds_spec_spectrum['albedo'], R=spectrum_resolution)
    else:
        w_noclouds, f_noclouds = jdi.mean_regrid(noclouds_spec_spectrum['wavenumber'],
                          noclouds_spec_spectrum['fpfs_reflected'], R=spectrum_resolution)
    
    
    ### Cloud-y spectrum:
    pdict, sdict, cdict, clouddict = pickle.load(open(savefiledirectory+'/cloudy-model-inputs.pkl','rb'))
    clouds_spec = pickle.load(open(savefiledirectory+'/cloudy-model.pkl','rb'))
    clouds_spec_spectrum = clouds_spec.spectrum(opa_mon, 
                    calculation='reflected', 
                    full_output=True)
    if plot_albedo:
        w_clouds, f_clouds = jdi.mean_regrid(clouds_spec_spectrum['wavenumber'],
                              clouds_spec_spectrum['albedo'], R=spectrum_resolution)
    else:
        w_clouds, f_clouds = jdi.mean_regrid(clouds_spec_spectrum['wavenumber'],
                          clouds_spec_spectrum['fpfs_reflected'], R=spectrum_resolution)
        
    pickle.dump([noclouds_spec_spectrum,1e4/w_noclouds, f_noclouds],
               open(savefiledirectory+'/cloud-free-spectrum-R'+str(spectrum_resolution)+'.pkl','wb'))
    pickle.dump([clouds_spec_spectrum,1e4/w_clouds, f_clouds],
               open(savefiledirectory+'/cloudy-spectrum-R'+str(spectrum_resolution)+'.pkl','wb'))
    
    # make plot:
    fig = plt.figure()
    plt.plot(1e4/w_clouds, f_clouds, color='black', label='Cloudy')
    plt.plot(1e4/w_noclouds, f_noclouds, color='darkcyan', label='Cloud-Free')
    plt.minorticks_on()
    plt.tick_params(axis='both',which='major',length =10, width=2,direction='in',labelsize=23)
    plt.tick_params(axis='both',which='minor',length =5, width=2,direction='in',labelsize=23)
    plt.xlabel(r"Wavelength [$\mu$m]", fontsize=25)
    plt.ylabel('Planet:Star Contrast', fontsize=25)
    plt.gca().set_yscale('log')
    plt.grid(ls=':')
    plt.legend(fontsize=15, loc='lower left')
    plt.tight_layout()
    plt.savefig('reflected-spectrum-plot.png', bbox_inches='tight')
    
    return noclouds_spec_spectrum, 1e4/w_noclouds, f_noclouds, clouds_spec_spectrum, 1e4/w_clouds, f_clouds, fig
    