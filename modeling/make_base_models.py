import numpy as np
from myastrotools.tools import MakeModelCloudFreePlanet, MakeModelCloudyPlanet, MakeModelCloudyAndCloudFreeSpectra
import astropy.units as u

## Planet:
planettype = 'ColdJup'
Tint = 100 # Internal Temperature of your Planet in K
Teq = 200 # planet equilibrium temperature 
radius = 1 #Rjup
msinij = 1
semi_major = 3

## Star:
T_star = 6000.0 # K, star effective temperature
logg = 4.437 #logg , cgs
metal = 0.0122 # metallicity of star
r_star = 0.8 # solar radius


## Climate:
nlevel = 91 # number of plane-parallel levels in your code
nofczns = 1 # number of convective zones initially. Let's not play with this for now.
nstr_upper = 85 # top most level of guessed convective zone
nstr_deep = nlevel -2 # this is always the case. Dont change this
nstr = np.array([0,nstr_upper,nstr_deep,0,0,0]) # initial guess of convective zones
rfacv = 0.5

## Opacities:
#
#planet_mh_str = '-100'# #log metallicity
planet_mh_str = '+050'#' #log metallicity
#planet_mh = float(planet_mh_str[1:])/100
planet_mh_CtoO_str = '100'#'0.5' # CtoO ratio


# save model params
directory = f'{planettype}-Teq{Teq}-Tint{Tint}-sep{semi_major}-mh{planet_mh_str}-co{planet_mh_CtoO_str}-Tstar{T_star}'
savefiledirectory = '/Volumes/Oy/Reflected-Light-Ames/models/'+directory

    
    
planet_properties = {
    'tint':Tint, 'Teq':Teq, 'radius':radius, 'radius_unit':u.Rjup,
     'mass':msinij, 'mass_unit': u.Mjup,
     'gravity': None, 'gravity_unit':None,
    'semi_major':semi_major, 'semi_major_unit': u.AU,
    'mh': planet_mh_str, 'CtoO':planet_mh_CtoO_str
}

star_properties = {
    'Teff':T_star, 'logg':logg, 'mh':metal, 'radius':r_star
}

climate_run_setup = {
    'nlevel':nlevel, 'nofczns':nofczns, 'nstr_upper':nstr_upper,
    'nstr_deep':nstr_deep, 'rfacv':rfacv
}


cj = MakeModelCloudFreePlanet(planet_properties, 
                        star_properties, 
                        use_guillotpt = True,
                        cdict = climate_run_setup,
                        climate_pbottom = 2,
                        climate_ptop = -6, 
                        savemodel = True,
                        savefiledirectory = savefiledirectory
             )

# cloud properties:
mean_molecular_weight = 2.2
fsed = 1
kz = 1e9
mh = 1 # Must always be 1 for virga runs
cloud_properties = {
    'fsed':fsed, 'mean_mol_weight': mean_molecular_weight, 'kz':kz, 'mh':mh
}
 
    
clouds, clouds_added = MakeModelCloudyPlanet(savefiledirectory, cloud_properties)

noclouds_spectrum, w_noclouds, f_noclouds, clouds_spectrum, w_clouds, f_clouds, fig = MakeModelCloudyAndCloudFreeSpectra(
                            savefiledirectory,
                            spectrum_wavelength_range = [0.5,1.8],
                            spectrum_calculation = 'reflected',
                            spectrum_resolution = 150,
                            calculation = "planet",
                            plot_albedo = False
                                      )

os.system('say "Done"')