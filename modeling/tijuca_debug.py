import picaso.justdoit as jdi
import pickle
local_ck_path = f'/Volumes/Oy/picaso/reference/kcoeff_2020/'
ck_db_name = local_ck_path + f'sonora_2020_feh+200_co_050_noTiOVO.data.196'
opacity_ck = jdi.opannection(ck_db=ck_db_name, wave_range = [0.4,2.0])

from virga import justdoit as vj
savefiledirectory = 'ReflectXGasGiantModelGrid-Neptune-Tstar3000-Rstar0.2-Teq187-sep0.1-rad0.35-mass0.05-mh100-co0.5-phase0'
pl, noclouds = pickle.load(open(savefiledirectory+'/cloud-free-model.pkl','rb'))

temperature = noclouds['temperature']
pressure = noclouds['pressure']
metallicity_TEMP = 0

# got molecules for cloud species:
mols = ['H2O','KCl', 'ZnS', 'NH3']
#mols = ['H2O','KCl', 'ZnS']
# add atmosphere from climate run:
atm = noclouds['ptchem_df']
# add kzz:
atm['kz'] = [1e9]*atm.shape[0]

clouds = jdi.inputs(calculation='planet') # start a calculation
pl.phase_angle(0)

# add gravity:
clouds.gravity(radius=0.35, radius_unit=u.Rjup, 
    mass = 0.05, mass_unit=u.Mjup)

clouds.star(opacity_ck, temp = 3000, metal = np.log10(1), logg = 5, 
    radius = 0.2, radius_unit=u.R_sun, 
    semi_major = 0.1, semi_major_unit = u.au, database='phoenix')
# add atmosphere from climate run with kzz:
clouds.atmosphere(df=atm)
# get clouds from reference:

meiff_directory ='/Volumes/Oy/virga/virga/reference/RefIndexFiles'

clouds_added = clouds.virga(mols, meiff_directory, fsed=0.1, mh=10**(metallicity_TEMP),
                     mmw = 2.2, full_output=True)

print('Done')
