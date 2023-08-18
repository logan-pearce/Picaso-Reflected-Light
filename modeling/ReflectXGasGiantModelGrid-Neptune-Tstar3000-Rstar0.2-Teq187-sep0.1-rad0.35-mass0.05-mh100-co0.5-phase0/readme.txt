
    ReflectX Gas Giant Picaso model 
    -------------------------------
    generated 2023-8-17-17:14:29 in python 3.8 
    Consult www.loganpearcescience/reflectx.html for more details 

    Model: 
     - Planet type: Neptune 
     - Star T_eff [K]: 3000 
     - Star radius [Rsun]: 0.2 
     - Planet Eq Temp [K]: 187 
     - Planet semi-major axis [au]: 0.1 
     - Planet radius [Rjup]: 0.35 
     - Planet mass [Mjup]: 0.05 
     - Planet metallicity: 100 
     - Planet C/O ratio: 0.5 
     - Phase: 0 

    Picaso input parameters: 
    These dictionaries are in machine-readable format (pickle) in '/cloud-free-model-inputs.pkl' 
    To read them in, run: `pdict, sdict, cdict = pickle.load(open(path_to_model+'/cloud-free-model-inputs.pkl','rb'))` 
     - Planet Properties Dictionary: {'tint': 100, 'Teq': 187.1302769665348, 'radius': 0.35, 'radius_unit': Unit("jupiterRad"), 'mass': 0.05, 'mass_unit': Unit("jupiterMass"), 'gravity': None, 'gravity_unit': None, 'semi_major': 0.1, 'semi_major_unit': Unit("AU"), 'mh': 100, 'CtoO': 0.5, 'phase': 0, 'num_tangle': 1, 'num_gangle': 6, 'noTiOVO': True, 'planet_mh_str': 200, 'local_ck_path': '/Volumes/Oy/picaso/reference/kcoeff_2020/'} 
     - Star Properties Dictionary: {'Teff': 3000, 'logg': 5.0, 'mh': 1, 'radius': 0.2, 'flux_unit': 'ergs cm^-2 s^-1 cm^-1'} 
     - Climate Run Parameters Dictionary: {'climate_pbottom': 2, 'climate_ptop': -6, 'nlevel': 91, 'nofczns': 1, 'nstr_upper': 85, 'nstr_deep': 89, 'rfacv': 0.5} 
     - Spectrum Calculation Dictionary: {'opacity_db': None, 'wave_range': [0.4, 2.0], 'calculation': 'reflected', 'R': 2000} 

    Virga-recommended condensation gases: ['H2O', 'KCl', 'ZnS'] 

    DIRECTORY CONTENTS:
    ------------------

    Cloud-free full PICASO output located in '/cloud-free-model.pkl'
    To read it in: `model = pickle.load(open(path_to_model+'/cloud-free-model.pkl','rb'))`

    Cloud-free spectrum located in '/cloud-free-spectrum-full-output-R'+str(specdict['R'])': 
     - Machine readable pickle file: `spectrum_dictionary = pickle.load(open(path_to_model+'/cloud-free-spectrum-full-output-R'+str(specdict['R']+'.pkl'), 'rb'))`
     - Human and machine readable csv: `spectrum_df = pandas.read_csv(path_to_model+'/cloud-free-spectrum-full-output-R'+str(specdict['R'])+'.csv', delim_whitespace=True)`

    `/terminal-output.txt`: A text file containing terminal output during the run. Useful for assessing
     model convergence.

    Plots:

     - `4SpectrumPlot`: plot of planet albedo spectrum, planet contrast (Fp/Fs), stellar spectrum, and planet spectrum (Fp/Fs * Fs)
     - `Abundances`: molecular abundances of the 8 most abundant molecules as a function of pressure (altitude)
     - `PTprofile`: Pressure/Temperature profile from climate calculation
     - `recommended-gases.html`: plot of PT profile against gas condensation curves with the gases used in the spectrum identified with heavy lines

    