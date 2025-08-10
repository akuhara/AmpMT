# AmpMT
## A modern fortran program to estimate normalized seismic moment tensor with non-double couple component from P-wave polarities and amplitudes. 

## Requirements
* GNU Fortran Compiler: Verified to work with version 12.3.0
* Opem MPI: Verified to work with version 4.1.2

## Compile 
* `cd src`
* `make`
It may rquire user-edit on the `src/Makefile` for appropriate paths to libararies etc.

## Usage 
`mpirun -np 20 mtamp [parameter file]`

## Parameter file example
```
polarity_file = obs.dat
station_file = /mnt/SunDisk/Noto_2024/Data/station.list
n_procs = 20
n_iter = 500000
n_burn = 250000
n_interval = 5000
n_chains = 5
n_cool = 1
temp_high = 300.0
dc_only = F
sample_prior = F
use_amp = T
```

## Polarity file example
```
75 240129.162657 37.89743 137.71090 11.303
  N_OBS  29
S18A 67.5 160.4 U 0.127270788244829
S14A 222.3 139.0 U 0.922355355570759
S13A 304.6 135.2 D -0.560374488663066
S15A 169.4 124.9 U 0.0488284060526949
S17A 19.3 122.4 U -1.0623461853513
S19A 100.4 120.4 D -0.0267965193199355
S21A 64.8 111.5 D -0.669528277096184
S10A 226.9 110.0 U -0.481746686896991
S24A 208.7 107.5 U -0.107461112492572
S16A 161.2 104.7 U -0.669806241046732
S20A 127.9 104.5 D 0.738411258712912
S22A 78.8 101.3 D -0.0271712813585943
S25A 193.7 102.6 U -0.448624526465005
S23A 104.5 100.5 D 0.787277221649404
S08A 206.6 98.4 U -0.557756356316013
S12A 183.1 97.0 U -1.11955683723866
S07A 228.1 94.7 D -0.957494168392587
S09A 196.3 88.4 U -1.0305967404544
S05A 213.6 87.8 U -1.01341928043107
S04A 228.8 84.6 D -0.44285086818451
S06A 204.0 67.9 U -0.657880232470654
S02A 217.4 67.9 U -0.77978223211326
S01A 229.1 67.8 U -0.883623095938827
S03A 209.0 67.8 U -0.776999293304033
E.YUOS 217.3 67.8 U -1.86012920653508
N.SUZH 222.5 67.8 U -1.33434140238179
S00A 234.0 67.8 U -0.927943393801674
SUZU 212.1 67.8 U -1.02354119223649
E.IDES 218.2 67.7 U -0.133095053723354
```
1st line: [ID]  [yymmdd.HHMMSS] [Lat.] [Lon.] [Dep.]

2nd line: [N_OBS] [# of stations]

3rd line or later: [Station name] [Azimuth] [Takeoff angle] [P-wave polarity (U or D)] [Log amplitude of P-wave]

NOTE: P-wave amplitude must be corrected for radiation pattern beforehand (need to remove local site effects, attenuation along propagation path etc.)

### Output files

## summary_00001.dat
Each line contains a model sampled by MCMC.

[Log-likelihood], [Mrr], [Mtt], [Mff], [Mrt], [Mrf], [Mtf], [T-strike], [T-dip], [B-strike], [B-dip], [P-strike], [P-dip], [h on ternary diagram], [v on ternary diagram], [h on source-type plot], [v on source-type plot], [normalization constant]



