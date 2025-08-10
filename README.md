# AmpMT
## A modern fortran program to estimate normalized seismic moment tensor with non-double couple component from P-wave polarities and amplitudes. 

## Requirements
* GNU Fortran Compiler: Verified to work with version 12.3.0
* Open MPI: Verified to work with version 4.1.2

## Compile 
* `cd src`
* `make`
Note: You may need to edit `src/Makefile` to specify appropriate paths to libraries and other dependencies.

## Usage 
`mpirun -np 20 mtamp [parameter file]`

## Example parameter file
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

## Example polarity file
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
File Format

1st line: [Event ID]  [yymmdd.HHMMSS] [Lat.] [Lon.] [Dep.]

2nd line: [N_OBS] [Number of stations]

3rd line and later: [Station name] [Azimuth] [Takeoff angle] [P-wave polarity (U/D)] [Log amplitude of P-wave]

NOTE: P-wave amplitudes must be corrected for the radiation pattern beforehand (i.e., remove local site effects, attenuation along the propagation path etc.)

### Example station file 
```
E.IDES 37.44005 137.25987 23
E.YUOS 37.51477 137.34449 26
N.SUZH 37.5266 137.2844 -154
S00A 37.599173 137.196407 -93
S01A 37.618877 137.307120 -103
S02A 37.588243 137.413705 -197
S03A 37.499645 137.434060 -111
S04A 37.676588 137.393687 -116
S05A 37.645897 137.500352 -365
S06A 37.557640 137.520065 -978
S07A 37.733918 137.481108 -238
S08A 37.703357 137.588205 -708
S09A 37.615397 137.606560 -1348
S10A 37.792092 137.568768 -418
S12A 37.672407 137.695132 -1140
S13A 37.938550 137.635842 -1597
S14A 37.849492 137.655698 -1275
S15A 37.799720 137.733750 -983
S16A 37.729863 137.782627 -1634
S17A 37.996615 137.755077 -1643
S18A 37.907917 137.743073 -1520
S19A 37.877072 137.849167 -1595
S20A 37.788047 137.887567 -1801
S21A 37.958532 137.875622 -1815
S22A 37.934910 137.953448 -1788
S23A 37.846537 137.956185 -1758
S24A 37.752425 137.610867 -617
S25A 37.711637 137.653563 -712
SUZU 37.4503 137.3587 10
```
File Format

[Station name] [Lat.] [Lon.] [Elevation (m)]


### Output files

## summary_00001.dat
Each line contains a model sampled by MCMC, with the following format:

[Log-likelihood], [Mrr], [Mtt], [Mff], [Mrt], [Mrf], [Mtf], [T-strike], [T-dip], [B-strike], [B-dip], [P-strike], [P-dip], [h on ternary diagram], [v on ternary diagram], [h on source-type plot], [v on source-type plot], [normalization constant]



