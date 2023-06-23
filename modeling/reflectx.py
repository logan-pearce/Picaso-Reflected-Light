import numpy as np
import astropy.units as u
import astropy.constants as c
import os
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt


def ConvertPlanetMHtoCKStr(m):
    prefixsign = np.sign(m)
    if prefixsign == -1:
        prefix = '-'
    else:
        prefix = '+'
    m = int(np.abs(m))
    if m / 100 >= 1.0:
        m = prefix+str(m)
    elif m == 0:
        m = '+000'
    else:
        m = prefix+'0'+str(m)
    return m

def ConvertCtoOtoStr(c):
    if c >= 1:
        cc = str(c*100).replace('.0','')
    else:
        cc = '0'+str(c*100).replace('.0','')
    return cc

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

def MakePTProflePlot(atm_df):
    #### Plot PT profile:
    fig = plt.figure(figsize=(10,10))
    plt.ylabel("Pressure [Bars]", fontsize=25)
    plt.xlabel('Temperature [K]', fontsize=25)
    plt.ylim(500,1e-6)
    plt.xlim(0,3000)
    plt.semilogy(atm_df['temperature'],
                atm_df['pressure'],color="r",linewidth=3)
    plt.minorticks_on() 
    plt.tick_params(axis='both',which='major',length =30, width=2,direction='in',labelsize=23)
    plt.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)
    #plt.legend(fontsize=15)

    return fig

def GetPlanetFlux(PlanetWNO,PlanetFPFS,StarWNO,StarFlux):
    from scipy.interpolate import interp1d
    func = interp1d(StarWNO,StarFlux,fill_value="extrapolate")
    ResampledStarFlux = func(PlanetWNO)
    return ResampledStarFlux*PlanetFPFS

def Make4SpectrumPlot(PlanetWNO,PlanetAlbedo,PlanetFPFS,StarWNO,StarFlux,colormap = 'magma',
                      plotstyle='magrathea'):
    import matplotlib
    plt.style.use(plotstyle)
    cmap = matplotlib.cm.get_cmap(colormap)
    n = 6
    cs = np.linspace(0,1,n)
    colors = cmap(cs)
    
    fig, axs = plt.subplots(2, 2, figsize=(12,10))
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[1, 0]
    ax4 = axs[1, 1]
    
    ax1.minorticks_on()
    ax1.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
    ax1.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)
    ax1.plot(1e4/PlanetWNO, PlanetAlbedo, lw=2, color=colors[1])
    ax1.set_yscale('log')
    #ax1.set_xlabel(r'Wavelength [$\mu$m]')
    ax1.set_ylabel('Planet Albedo')
    ax1.grid(ls=':')
    
    ax2.minorticks_on()
    ax2.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
    ax2.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)
    ax2.plot(1e4/PlanetWNO, PlanetFPFS, lw=2, color=colors[2])
    ax2.set_yscale('log')
    #ax2.set_xlabel(r'Wavelength [$\mu$m]')
    ax2.set_ylabel('Planet:Star contrast')
    ax2.grid(ls=':')
    
    ax3.minorticks_on()
    ax3.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
    ax3.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)
    ax3.plot(1e4/StarWNO, StarFlux, lw=2, color=colors[3])
    ax3.set_yscale('log')
    ax3.set_xlabel(r'Wavelength [$\mu$m]')
    ax3.set_ylabel(r'Star flux [ergs cm$^{-2}$ s$^{-1}$ cm$^{-1}$]')
    ax3.grid(ls=':')
    
    PlanetFlux = GetPlanetFlux(PlanetWNO,PlanetFPFS,StarWNO,StarFlux)
    ax4.minorticks_on()
    ax4.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
    ax4.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)
    ax4.plot(1e4/PlanetWNO, PlanetFlux, lw=2, color=colors[4])
    ax4.set_yscale('log')
    ax4.set_xlabel(r'Wavelength [$\mu$m]')
    ax4.set_ylabel(r'Planet flux [ergs cm$^{-2}$ s$^{-1}$ cm$^{-1}$]')
    ax4.grid(ls=':')
    plt.tight_layout()
    return fig

def ComputeSpectrum(atm_df, pdict, sdict, specdict, calculation = 'planet', 
                   clouddict = None,
                   saveplots = True, savefiledirectory = None, 
                   ComputePlanetFlux = True, make4spectrumplot = True,
                   MixingRatioPlot = True, AbundancesPlot = True, n_mols_to_plot = 8,
                   PTplot = True,
                   plotstyle = 'magrathea'):
    
    import picaso.justdoit as jdi
    opacity_db = specdict['opacity_db']
    if opacity_db == None:
        opa_mon = jdi.opannection(wave_range=specdict['wave_range'])
    else:
        opa_mon = jdi.opannection(filename_db = opacity_db, wave_range=specdict['wave_range'])
    spec = jdi.inputs(calculation=calculation)

    spec.phase_angle(phase=pdict['phase']*np.pi/180, num_tangle=pdict['num_tangle'],
                     num_gangle=pdict['num_gangle'])
    if not pdict['gravity']:
        spec.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        spec.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])

    # set up star:
    spec.star(opa_mon, temp = sdict['Teff'], metal = np.log10(sdict['mh']), logg = sdict['logg'], 
            radius = sdict['radius'], radius_unit = u.R_sun, 
            semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'], database = 'phoenix')

    spec.atmosphere(df=atm_df)

    if clouddict != None:
        clouds_added = spec.virga(
            clouddict['condensates'], 
            clouddict['virga_mieff'], 
            fsed=clouddict['fsed'], mh=10**(0),
            mmw = clouddict['mean_mol_weight'], full_output=True)

    spec_df = spec.spectrum(opa_mon, 
                            calculation='reflected', 
                            full_output=True)
    wno, alb, fpfs, full_output = spec_df['wavenumber'], spec_df['albedo'], \
                                    spec_df['fpfs_reflected'],  spec_df['full_output']
    wno,fpfs = jdi.mean_regrid(spec_df['wavenumber'],
                          spec_df['fpfs_reflected'], R=specdict['R'])
    wno,alb = jdi.mean_regrid(spec_df['wavenumber'],
                          spec_df['albedo'], R=specdict['R'])
    
    if '/' not in savefiledirectory:
        savefiledirectory = savefiledirectory+'/'

    if saveplots:
        plt.style.use(plotstyle)
        # Albedo plot
        plt.figure(figsize=(8,6))
        plt.minorticks_on()
        plt.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
        plt.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)

        plt.plot(1e4/wno, alb, lw=2)

        plt.gca().set_xlabel(r'Wavelength [$\mu$m]')
        plt.gca().set_ylabel('Geo Albedo')
        plt.grid(ls=':')
        plt.tight_layout()
        plt.savefig(savefiledirectory+'-albedo-spectrum-R'+str(specdict['R'])+'.png',
                    bbox_inches='tight')

        # reflected plot
        plt.figure(figsize=(8,6))
        plt.minorticks_on()
        plt.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
        plt.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)

        plt.plot(1e4/wno, fpfs, lw=2)
        plt.gca().set_yscale('log')

        plt.gca().set_xlabel(r'Wavelength [$\mu$m]')
        plt.gca().set_ylabel('Planet:Star contrast')
        plt.grid(ls=':')
        plt.tight_layout()
        plt.savefig(savefiledirectory+'-fpfs-spectrum-R'+str(specdict['R'])+'.png',
                    bbox_inches='tight')
        
    if make4spectrumplot:
        StarWNO,StarFlux = spec.inputs['star']['wno'],spec.inputs['star']['flux']
        fig = Make4SpectrumPlot(wno,alb,fpfs,StarWNO,StarFlux,colormap = 'magma')
        fig.savefig(savefiledirectory+'-4SpectrumPlot-R'+str(specdict['R'])+'.png',
                    bbox_inches='tight')
        
    if MixingRatioPlot:
        import picaso.justplotit as jpi
        mrp = jpi.mixing_ratio(full_output, plot_height=900, plot_width=600)
        # Save mixing ratio plot for future inspection
        from bokeh.plotting import output_file, save
        output_file(savefiledirectory+"-mixing-ratios.html")
        save(mrp)

    if AbundancesPlot:
        sort = np.argsort(atm_df.loc[0])[::-1]
        highest_abundance = np.array(atm_df.keys()[sort])
        highest_abundance = np.delete(highest_abundance,np.where(highest_abundance=='temperature'))
        highest_abundance = np.delete(highest_abundance,np.where(highest_abundance=='pressure'))
        highest_abundance = np.delete(highest_abundance,np.where(highest_abundance=='e-'))
        import matplotlib
        cmap = matplotlib.cm.get_cmap('Spectral')
        mols = highest_abundance[:n_mols_to_plot+1]
        linestyles = ['-','-.']*len(mols)
        lineweights = [4,4,2,2]*len(mols)
        n = len(mols)
        cs = np.linspace(0.1,1,n)
        colors = cmap(cs)
        plt.figure(figsize=(10,10))
        plt.minorticks_on()
        plt.tick_params(axis='both',which='major',length =30, width=2,direction='in',labelsize=23)
        plt.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)

        plt.ylabel("Pressure [Bars]", fontsize=25)
        plt.xlabel('Abundance', fontsize=25)
        plt.ylim(500,1e-4)
        plt.xlim(left=1e-6)

        for i,mol in enumerate(mols):
            plt.plot(atm_df[mol],
                    atm_df['pressure'],
                    color=colors[i],linewidth=lineweights[i],label=mol,ls=linestyles[i])
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        #plt.gca().invert_yaxis()
        plt.legend(fontsize=20)
        plt.savefig(savefiledirectory+'-abundances.png',
                    bbox_inches='tight')
        
    if PTplot:
        fig = MakePTProflePlot(atm_df)
        fig.savefig(savefiledirectory+'-PTprofile.png',
                    bbox_inches='tight')
        
    if ComputePlanetFlux:
        PlanetFlux = GetPlanetFlux(wno, fpfs, StarWNO, StarFlux)
        # reflected plot
        plt.figure(figsize=(8,6))
        plt.minorticks_on()
        plt.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
        plt.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)

        plt.plot(1e4/wno, PlanetFlux, lw=2)
        plt.gca().set_yscale('log')

        plt.gca().set_xlabel(r'Wavelength [$\mu$m]')
        plt.gca().set_ylabel(r'Planet flux [ergs cm$^{-2}$ s$^{-1}$ cm$^{-1}$]')
        plt.grid(ls=':')
        plt.tight_layout()
        plt.savefig(savefiledirectory+'-planet-spectrum-R'+str(specdict['R'])+'.png',
                    bbox_inches='tight')
        return wno, alb, fpfs, full_output, PlanetFlux, StarWNO, StarFlux

    return wno, alb, fpfs, full_output


def MakeModelCloudFreePlanet(pdict, sdict,
                calculation = "planet",
                use_guillotpt = True,
                user_supplied_ptprofile = None,
                cdict = None,
                compute_spectrum = True,
                specdict = None,
                savefiledirectory = None
             ):
    
    ''' Wrapper for PICASO functions for building a cloud-free planet base model
    Args:
        pdict (dict): dictionary of planet parameter inputs
        sdict (dict): dictionary of star parameter inputs
        calculation (str): picaso input for object, "planet" or "brown dwarf"
        use_guillotpt (bool): if True, use Guillot PT approximation. Else user must supply initial PT profile
        user_supplied_ptprofile (df): user supplied pt profile for picaso
        cdict (dict): dictionary of climate run setup params
        specdict (dict): dictionary of spectrum setup params
        savefilename (str): filename and path for the model to be saved.

    Returns:
        pl: picaso planet model inputs
        noclouds: picaso object after climate run before clouds

    Writes to disk:
        terminal_output.txt: text file containing terminal output from run
        recommended-gases.html: bokeh plot fro virga of recommended molecules for cloud run
        cloud-free-model.pkl: pickle file of jdi inputs object (pl) and atmosphere (noclouds)
        cloud-free-model-inputs.pkl: pickle file of dictionary inputs pdict, sdict, cdict
        cloud-free-spectrum-full-output-RXXX.pkl: pickle file of dictionary of all spectrum objects - 
                wavenumber, albedo spectrum, planet/star flux contrast, planet flux, star spectrum wavenumber, star flux
        cloud-free-fpfs-spectrum-RXXX.png: plot of planet star flux contrast
        cloud-free-alb-spectrum-RXXX.png: plot of planet albedo spectrum
        cloud-free-4SpectrumPlot-RXXX.png: plot of planet albedo, contrast, star flux, and planet flux
    '''
    import warnings
    warnings.filterwarnings('ignore')
    import picaso.justdoit as jdi
    
    import sys
    import os
    # Make directory to store run results:
    os.system('mkdir '+savefiledirectory)
    
    # Set all terminal output to write to file:
    f = open(savefiledirectory+"/terminal_output.txt", 'w')
    sys.stdout = f
    # Additional info for xarray (not used)
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
    PlanetCO = ConvertCtoOtoStr(pdict['CtoO'])
    PlanetMHStr = ConvertPlanetMHtoCKStr(pdict['planet_mh_str']).replace('.0','')

    if pdict['noTiOVO']:
        ck_db_name = pdict['local_ck_path'] + f'sonora_2020_feh{PlanetMHStr}_co_{PlanetCO}_noTiOVO.data.196'
    else:
        ck_db_name = pdict['local_ck_path'] + f'sonora_2020_feh{PlanetMHStr}_co_{PlanetCO}.data.196'
    print(ck_db_name)
    # Set opacity connection:
    opacity_ck = jdi.opannection(ck_db=ck_db_name, wave_range = specdict['wave_range'])
    
    # initialize model:
    pl = jdi.inputs(calculation= calculation, climate = True)
    
    ### set up planet:
    # input effective temperature
    pl.effective_temp(pdict['tint']) 
    # add gravity:
    if not pdict['gravity']:
        pl.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        pl.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
        
    # set up star:
    pl.star(opacity_ck, temp = sdict['Teff'], metal = np.log10(sdict['mh']), logg = sdict['logg'], 
            radius = sdict['radius'], radius_unit = u.R_sun, 
            semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'], database = 'phoenix')
    
    # add phase:
    phase = pdict['phase']
    # If full phase:
    if phase == 0:
        # Use symmetry to speed up calculation.
        num_tangle = 1
        pdict['num_tangle'] = 1
    else:
        num_tangle = pdict['num_tangle']
    pl.phase_angle(phase=phase*np.pi/180, num_tangle=num_tangle, num_gangle=pdict['num_gangle'])
    
    # Add initial P(T) profile guess:
    if use_guillotpt:
        pt = pl.guillot_pt(pdict['Teq'], nlevel=cdict['nlevel'], T_int = pdict['tint'], 
                              p_bottom=cdict['climate_pbottom'], p_top=cdict['climate_ptop'])
    else:
        pt = user_supplied_ptprofile

    # initial PT profile guess:
    temp_guess = pt['temperature'].values 
    press_guess = pt['pressure'].values
    # Input climate params:
    nstr = np.array([0,cdict['nstr_upper'],cdict['nstr_deep'],0,0,0]) # initial guess of convective zones
    # Set up climate run inputs:
    pl.inputs_climate(temp_guess= temp_guess, pressure= press_guess, 
                  nstr = nstr, nofczns = cdict['nofczns'] , rfacv = cdict['rfacv'])
    print('starting climate run')
    # Compute climate:
    noclouds = pl.climate(opacity_ck, save_all_profiles=True, with_spec=True)
    # Set atm to climate run results:
    pl.atmosphere(df=noclouds['ptchem_df'])

    # Generate recommended molecules for clouds:
    from virga import justdoit as vj
    # pressure temp profile from climate run:
    temperature = noclouds['temperature']
    pressure = noclouds['pressure']
    # metallicity = pdict['mh'] #atmospheric metallicity relative to Solar
    metallicity_TEMP = 0 # mh must be 1 for this part
    mmw = 2.2
    # got molecules for cloud species:
    recommended, fig = vj.recommend_gas(pressure, temperature, 10**(metallicity_TEMP), 
                                   mmw, plot=True, outputplot = True)
    print('Virga recommended molecules:', recommended)
    import pickle
    pickle.dump(recommended, open(savefiledirectory+'/virga-recommended-molecules.pkl','wb'))
    
    # Save recommended plot for future inspection
    from bokeh.plotting import figure, output_file, save
    output_file(savefiledirectory+"/recomended-gasses.html")
    save(fig)

    sdict.update({'flux_unit':'ergs cm^-2 s^-1 cm^-1'})

    # Pickle info and save:
    pickle.dump([pl, noclouds], open(savefiledirectory+'/cloud-free-model.pkl','wb'))
    pickle.dump([pdict, sdict, cdict], open(savefiledirectory+'/cloud-free-model-inputs.pkl','wb'))

    # If you want to compute spectrum here (recommended)
    if compute_spectrum:
        ### Make cloud-free spectrum:
        wno, alb, fpfs, full_output, PlanetFlux, StarWNO, StarFlux = ComputeSpectrum(noclouds['ptchem_df'], 
                                                                  pdict, sdict, specdict, 
                                                                  calculation = 'planet',
                                                                  saveplots = True, 
                                                                  clouddict = None,
                                                                  savefiledirectory = savefiledirectory+'/cloud-free', 
                                                                  ComputePlanetFlux = True, 
                                                                  make4spectrumplot = True,
                                                                  AbundancesPlot = True, 
                                                                  MixingRatioPlot = True)
        # dump it out:
        outdict = {'wavenumber':wno, 'albedo spectrum':alb,
                   'planet star contrast':fpfs, 'full_output':full_output,
                   'planet flux': PlanetFlux, 'star wavenumber':StarWNO, 'star flux':StarFlux}
        pickle.dump(outdict, 
                    open(savefiledirectory+'/cloud-free-spectrum-full-output-R'+str(specdict['R'])+'.pkl','wb'))
        from scipy.interpolate import interp1d
        func = interp1d(StarWNO,StarFlux,fill_value="extrapolate")
        ResampledStarFlux = func(wno)
        outdf = pd.DataFrame(data={'wavelength [um]':1e4/wno,
                                   'planet flux [ergs/cm2/s/cm]':PlanetFlux,
                                   'star flux [ergs/cm2/s/cm]':ResampledStarFlux,
                                   'albedo':alb,
                                   'planet-star contrast':fpfs }
                             )
        outdf.to_csv(savefiledirectory+'/cloud-free-spectrum-R'+str(specdict['R'])+'.csv', index=False, sep=' ')
        
    f.close()
    
    return pl, noclouds
    

def MakeModelCloudyPlanet(savefiledirectory, clouddict, 
                          specdict,
                          cloud_filename_prefix,
                          compute_spectrum = True,
                          saveplots = True,
                          calculation = 'planet'):
    
    ''' Wrapper for PICASO functions for building a planet model
    Args:
        savefiledirectory (str): directory containing picaso cloud-free model base case.
        clouddict (dict): dictionary of cloud parameter inputs
        cloud_filename_prefix (str): name to prepend output files
        calculation (str): picaso input for object, "planet" or "brown dwarf"

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
    
    #opa_mon = jdi.opannection()
    PlanetMHStr = ConvertPlanetMHtoCKStr(pdict['planet_mh_str'])
    PlanetCO = ConvertCtoOtoStr(pdict['CtoO'])
    if pdict['noTiOVO']:
        ck_db_name = pdict['local_ck_path'] + f'sonora_2020_feh{PlanetMHStr}_co_{PlanetCO}_noTiOVO.data.196'
    else:
        ck_db_name = pdict['local_ck_path'] + f'sonora_2020_feh{PlanetMHStr}_co_{PlanetCO}.data.196'
    print(ck_db_name)
    # Set opacity connection:
    opacity_ck = jdi.opannection(ck_db=ck_db_name, wave_range = specdict['wave_range'])

    from virga import justdoit as vj
    # pressure temp profile from climate run:
    temperature = noclouds['temperature']
    pressure = noclouds['pressure']
    #metallicity = pdict['mh'] #atmospheric metallicity relative to Solar
    metallicity_TEMP = 0
    # got molecules for cloud species:
    if clouddict['condensates'] == 'virga recommend':
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
        mols = clouddict['condensates']
    # add atmosphere from climate run:
    atm = noclouds['ptchem_df']
    # add kzz:
    atm['kz'] = [clouddict['kz']]*atm.shape[0]

    # Just set all this up again:
    clouds = jdi.inputs(calculation=calculation) # start a calculation
    
    # add phase:
    phase = pdict['phase']
    # If full phase:
    if phase == 0:
        # Use symmetry to speed up calculation.
        num_tangle = 1
        pdict['num_tangle'] = 1
    else:
        num_tangle = pdict['num_tangle']
    pl.phase_angle(phase=phase*np.pi/180, num_tangle=num_tangle, num_gangle=pdict['num_gangle'])

    # add gravity:
    if not pdict['gravity']:
        clouds.gravity(radius=pdict['radius'], radius_unit=pdict['radius_unit'], 
            mass = pdict['mass'], mass_unit=pdict['mass_unit'])
    else:
        clouds.gravity(gravity=pdict['gravity'], gravity_unit=pdict['gravity_unit'])
    # add star:
    # clouds.star(opa_mon, temp = sdict['Teff'], metal = np.log10(sdict['mh']), logg = sdict['logg'], 
        # radius = sdict['radius'], radius_unit=u.R_sun, 
        # semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'], database='phoenix')
    clouds.star(opacity_ck, temp = sdict['Teff'], metal = np.log10(sdict['mh']), logg = sdict['logg'], 
        radius = sdict['radius'], radius_unit=u.R_sun, 
        semi_major = pdict['semi_major'], semi_major_unit = pdict['semi_major_unit'], database='phoenix')
    # add atmosphere from climate run with kzz:
    clouds.atmosphere(df=atm)
    # get clouds from reference:

    meiff_directory ='/Volumes/Oy/virga/virga/reference/RefIndexFiles'
    clouddict['virga_mieff'] = meiff_directory
    #meiff_directory = clouddict['virga_mieff']

    metallicity_TEMP = 0
    clouds_added = clouds.virga(mols, meiff_directory, fsed=clouddict['fsed'], mh=10**(metallicity_TEMP),
                         mmw = clouddict['mean_mol_weight'], full_output=True)

    clouddict.update({'condensates':mols})

    savefiledirectory = savefiledirectory+'/'
    outdir = savefiledirectory+cloud_filename_prefix+'/'
    os.system('mkdir '+ outdir)

    pickle.dump([clouds, clouds_added],
                    open(outdir+cloud_filename_prefix+'-cloudy-model.pkl','wb'))
    pickle.dump([pdict, sdict, cdict, clouddict],
                open(outdir+cloud_filename_prefix+'-cloudy-model-inputs.pkl','wb'))
    
    # If you want to compute spectrum here (recommended)
    if compute_spectrum:
        ### Make cloud-free spectrum:
        wno, alb, fpfs, full_output, PlanetFlux, StarWNO, StarFlux = ComputeSpectrum(noclouds['ptchem_df'], 
                                                                  pdict, sdict, specdict, 
                                                                  clouddict = clouddict,
                                                                  calculation = 'planet',
                                                                  saveplots = True, 
                                                                  savefiledirectory = outdir+cloud_filename_prefix, 
                                                                  ComputePlanetFlux = True, 
                                                                  make4spectrumplot = True,
                                                                  MixingRatioPlot = False,
                                                                  AbundancesPlot = False,
                                                                  PTplot = False)
        # dump it out:
        outdict = {'wavenumber':wno, 'albedo spectrum':alb,
                   'planet star contrast':fpfs, 'full output':full_output,
                   'planet flux': PlanetFlux, 'star wavenumber':StarWNO, 'star flux':StarFlux}
        pickle.dump(outdict, 
                    open(outdir+cloud_filename_prefix+'-spectrum-full-output-R'+str(specdict['R'])+'.pkl','wb'))
        from scipy.interpolate import interp1d
        func = interp1d(StarWNO,StarFlux,fill_value="extrapolate")
        ResampledStarFlux = func(wno)
        outdf = pd.DataFrame(data={'wavelength [um]':1e4/wno,
                                   'planet flux [ergs/cm2/s/cm]':PlanetFlux,
                                   'star flux [ergs/cm2/s/cm]':ResampledStarFlux,
                                   'albedo':alb,
                                   'planet-star contrast':fpfs }
                             )
        outdf.to_csv(outdir+cloud_filename_prefix+'-spectrum-R'+str(specdict['R'])+'.csv', index=False, sep=' ')
    f.close()

    return clouds, clouds_added



def RegridSpectrum(directory, spectrum_full_output_pickle_file, R, make4spectrumplot = False):
    ''' Regrid a spectrum
    
    Args:
        directory (str): location of filder containing orginal FULL OUTPUT spectrum pickle file
        spectrum_full_output_pickle_file (str): pickle file of dataframe of original spectrum. 
                    Ex: /cloud-free-spectrum-full-output-R2000.pkl
        R (flt): R of new spectrum
    
    '''
    import pickle
    import os
    import picaso.justdoit as jdi

    df = pickle.load(open(directory+spectrum_full_output_pickle_file,'rb'))
    wno, alb, fpfs, full_output, PlanetFlux, StarWNO, StarFlux = df['wavenumber'], df['albedo spectrum'], df['planet star contrast'],\
        df['full_output'], df['planet flux'], df['star wavenumber'], df['star flux']
    
    if 'cloud-free' in spectrum_full_output_pickle_file:
        fileout = 'cloud-free-'
    elif 'cloudy' in spectrum_full_output_pickle_file:
        fileout = 'cloudy-'

    if make4spectrumplot:
        fig = Make4SpectrumPlot(wno,alb,fpfs,StarWNO,StarFlux,colormap = 'magma')
        fig.savefig(directory+fileout+'-4SpectrumPlot-R'+str(R)+'.png',
                    bbox_inches='tight')
    
    # Albedo plot
    Rwno, Ralb = jdi.mean_regrid(wno, alb , R = R)
    plt.figure()
    plt.figure(figsize=(8,6))
    plt.minorticks_on()
    plt.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
    plt.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)

    plt.plot(1e4/Rwno, Ralb, lw=2)

    plt.gca().set_xlabel(r'Wavelength [$\mu$m]')
    plt.gca().set_ylabel('Geo Albedo')
    plt.grid(ls=':')
    plt.tight_layout()
    plt.savefig(directory+fileout + 'albedo-spectrum-R' + R + '.png',bbox_inches='tight')

    # reflected plot
    Rwno, Rfpfs = jdi.mean_regrid(wno, fpfs , R = R)
    plt.figure()
    plt.figure(figsize=(8,6))
    plt.minorticks_on()
    plt.tick_params(axis='both',which='major',length =20, width=2,direction='in',labelsize=23)
    plt.tick_params(axis='both',which='minor',length =10, width=2,direction='in',labelsize=23)

    plt.plot(Rwno, Rfpfs, lw=2)
    plt.gca().set_yscale('log')

    plt.gca().set_xlabel(r'Wavelength [$\mu$m]')
    plt.gca().set_ylabel('Planet:Star contrast')
    plt.grid(ls=':')
    plt.tight_layout()
    plt.savefig(directory+fileout + 'fpfs-spectrum-R' + R + '.png',bbox_inches='tight')

    pickle.dump([wno, alb, fpfs, full_output], open(directory+fileout+'-spectrum-full-output-R' + str(R) + '.pkl','wb'))
    
    from scipy.interpolate import interp1d

    from scipy.interpolate import interp1d
    func = interp1d(StarWNO,StarFlux,fill_value="extrapolate")
    ResampledStarFlux = func(Rwno)
    RPlanetFlux = ResampledStarFlux*Rfpfs

    outdf = pd.DataFrame(data={'wavelength [um]':1e4/Rwno,
                                'planet flux [ergs/cm2/s/cm]':RPlanetFlux,
                                'star flux [ergs/cm2/s/cm]':ResampledStarFlux,
                                'albedo':Ralb,
                                'planet-star contrast':Rfpfs }
                            )
    outdf.to_csv(directory+fileout+cloud_filename_prefix+'-spectrum-R'+str(R)+'.csv', index=False, sep=' ')
    