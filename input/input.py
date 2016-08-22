ppw = 128 # points per wavelength

position = ppw/16
thickness = ppw/16
foil = 0.7

PI = 3.14159265358979

## MISCELLANEOUS

def Block(x, xmin, xmax):
	if x > xmin and x < xmax:
		return 1
	else:
		return 0

## GENERAL PARAMETERS

MAX_Z = 6*ppw # number of the steps for coordinate
MAX_T = 10*ppw+1 # number of the steps for time
dt = 2*PI/ppw # step size for time
dz = 2*PI/ppw # step size for coordinate
THETA = 0 #PI/4 # incident angle for laser pulse

## PLASMA PARAMETERS

N_0 = 30 # overcritical parameter
NUM_SP = 3 # number of species

## SPECIE 0
MASS_0   =   1 # mass of the specie (relative to electron one)
CHARGE_0 =   1 # charge of the specie (relative to electron one, i.e. electron charge equals to +1)
MAX_P_0  = 20000 # number of the steps for momentum
dp_0 = 0.01 # step size for momentum (relative to Mc)
T_init_0 = 0.1/511 # initial temperature (relative to rest energy)
MEAN_P_0 = 0. # mean momentum (in Mc)
def PROFILE_0 (z):
	return Block(z, 2*ppw, (2+foil)*ppw)

## SPECIE 1
MASS_1   = 55.8*1835.3 # mass of the specie (relative to electron one)
CHARGE_1 = -24 # charge of the specie (relative to electron one, i.e. electron charge equals to +1)
MAX_P_1  = 2000 # number of the steps for momentum
dp_1 = 0.002 # step size for momentum (relative to Mc)
T_init_1 = 2e-5 # initial temperature (relative to rest energy)
MEAN_P_1 = 0. # mean momentum (in Mc)
def PROFILE_1 (z):
	return Block(z, 2*ppw+position, 2*ppw+position+thickness)

## SPECIE 2
MASS_2   = 197*1835.3 # mass of the specie (relative to electron one)
CHARGE_2 = -68 # charge of the specie (relative to electron one, i.e. electron charge equals to +1)
MAX_P_2  = 2000 # number of the steps for momentum
dp_2 = 0.002 # step size for momentum (relative to Mc)
T_init_2 = 2e-5 # initial temperature (relative to rest energy)
MEAN_P_2 = 0. # mean momentum (in Mc)
def PROFILE_2 (z):
	return 0 #PROFILE_0 (z) - PROFILE_1 (z)

## FIXED IONS
def FIXED_IONS_PROFILE (z):
	return PROFILE_0 (z) - PROFILE_1 (z) # profile of the fixed ions concentration distribution

## PULSE PARAMETERS

A = 50 # electromagnetic pulse amplitude
t_rise = dt # pulse rising time (only for trapezeidal and sin squared slopes)
t_length = 5*2*PI # pulse duration
t_delay = t_length # delay before pulse appearing
source = 10 # the position of a source of an electromagnetic wave in respect to the PML-layer

def PULSE_X(t):
	#return A*Block(t, t_delay, t_delay + t_length)*sin(t - t_delay) # a slope of the laser pulse
	return A*exp(-sqr((t-t_delay)/t_length*2))*cos(t) # a slope of a laser pulse
def PULSE_Y(t):
	return A*exp(-sqr((t-t_delay)/t_length*2))*sin(t) # a slope of a laser pulse

## TEST PARTICLE PARAMETERS

NUM_PRT = 0
MASS_PRT = MASS_1
CHARGE_PRT = CHARGE_1
start_point = 3*ppw
interval = 1

## OUTPUT PARAMETERS

output_directory_name = str(foil)+'_'+str(position)+'_'+str(thickness) # name of an output directory

save_format = 'gzip' # format of saving files (txt, bin, gzip)
save_dt = ppw/16 # time interval of saving

save_fields = 1 # if fields saving in files is needed
save_fields_format = '' # format of saving files
save_fields_dt = 0 # time interval of saving

save_concs = 1 # if concentrations saving in files is needed
save_concs_format = '' # format of saving files
save_concs_dt = 0 # time interval of saving

save_dstr = 1 # if distribution functions storing in file is needed
save_dstr_format = '' # format of saving files
save_dstr_dt = MAX_T/2 # time interval of saving
