

## In this study, we use a hybrid resolution approach. We use $90\times52$ 

MITgcm configuration details (can be also found in *data*) 

**Gridding parameters**

 &PARM04
 usingSphericalPolarGrid = .TRUE.,
 delZ=10.,10.,15.,20.,20.,25.,35.,50.,75.,100.,150.,200.,275.,350.,415.,450.,500.,500.,500.,500.,500.,500.,500.,
 ygOrigin=-88.,
 delY=4*2.,19*4.,8*1.,21*4.,
 delX=90*4.,
 &

**Input datasets**

 &PARM05\
 bathyFile       = 'bathy_90x52_mod_f64.bin',\
 hydrogThetaFile = 'woa2018_tempr8fill90x52x23.bin',\
 hydrogSaltFile  = 'woa2018_saltr8fill90x52x23.bin',
 
### +++++++++++++++++++++++++++++++++++++++++++++++++ 

The important changes are done in the *radmod* parameters in PlaSim (puma_namelist) 

&RADPAR\
th2oc = 0.02\
tswr1   = 0.08\
tswr2   = 0.08\
tswr3   = 0.0055\
acllwr  = 0.05\
NDCYCLE =     1\
CO2     = 360.0\
DAWN    =     0.0000\
GSOL0   = 1367.0000

### ----------------------------------------------------------------------

Experiment steps

1. We first ran the model with the modified parameters stated above for 50 years from the climatological initial conditions. We call it "control (ctrl)" run. The output files are kept in EXP/CTRL along with the modified scripts and input files. 
