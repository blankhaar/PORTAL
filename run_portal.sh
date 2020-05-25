#!/usr/bin/env bash
#SBATCH -A C3SE508-18-3 -p astro1
#SBATCH -n 32
#SBATCH -t 4-00:00:00

###############################################################################
###############################################################################
#############################INPUT PARAMETERS##################################
###############################################################################
###############################################################################

#input file names

#input from LIME for simulation
pop='/priv/c3-astro-data1/PORTAL/pol_files/lime/agb/agb_12CO-v2_small.pop'
grid='/priv/c3-astro-data1/PORTAL/pol_files/lime/agb/agb_12CO-v2_small.vtk'
dust='/priv/c3-astro-data1/PORTAL/pol_files/jena_thin_e6.tab'
mol='/priv/c3-astro-data1/PORTAL//pol_files/CO_coro_v2.dat'

#make sure to not run two scripts with the same
#input file names in the same map
input_file_1='pol.dat'
input_file_2='finish_pol.dat'

#output map
output_map='output'

#magnetic field character
b_rad='0'	#radial character 
b_tor='1'	#toroidal character
b_pol='0'	#poloidal character
b_dip='1.d0'    #if not zero, then B=B_dip

#you want to run program in star mode?
star_mode=1	#=1 then in star mode
T_star='2.5E3'
R_star='2.99196E11'	#in meters
X_star='0.d0'	#x coordinate in meters
Y_star='0.d0'	#y coordinate in meters
Z_star='0.d0'	#z coordinate in meters

#dust to gass mass ratio
dust_ii_gas='0.01d0'

#integration parameters
nth=64
nph=64
nvel=10

#parameter-length parameters
n_depth=10000
n_depth_rt=10000
max_ne_temp=10000

#maximum irreducible tensor element (>6 is advised)
kmax=6

#NEW:
#be sure to define this one correctly ---otherwise bogus angular scale
#how far away is the source (in meters)?
s_d='3.086e16'

#extra info?
extra_info=1	#print out rates and J2/J0 per node

#ray-tracing info
#how many transitions do you want to ray trace
ntrans=3
#what are their assoicated levels
l1[0]=1
l1[1]=2
l1[2]=3

l2[0]=2
l2[1]=3
l2[2]=4

freqs[0]='115e9'
freqs[1]='230e9'
freqs[2]='345e9'

#how many global angles are going to be ray traced
ninc=2
#what are their associated angles (in radians)
th_inc[0]='0.d0'
th_inc[1]='0.78539816339d0'

ph_inc[0]='0.d0'
ph_inc[1]='0.d0'

#number of velocity channels (=2*nvel+1) and resolution and centre
nvel_raytrace=0
dvel_raytrace='50.e0'	#m/s
v0='0.d0'

#number of pixels (=2*im_pix+1) and size of image (relative to simulation)
pix_mod=1       #if 1, then we'll do a pixel ray-tracing 
im_pix=100
im_size='0.7'

#NEW: this one not needed in new version
delta_pix="2.7778e-4"

#NEW: define fits
fits=1    #if 1 ---> fits file output


###############################################################################
###############################################################################
#################################JOB SCRIPT####################################
###############################################################################
###############################################################################

if [ ! -d "$output_map" ]; then
  mkdir $output_map
fi

cat << EOF > $input_file_1

$star_mode $T_star $R_star $X_star $Y_star $Z_star
$nth $nph $nvel
$dust_ii_gas
$b_rad $b_tor $b_pol $b_dip
$kmax $extra_info 
$max_ne_temp $n_depth

EOF

#input file for the ray-tracing gets an additional line where we have the 
#number 'fits' and the number 's_d'. 
#fits=1 ---> output files in fits
#s_d	---> source distance in meters (used for the conversion to angular scale)

cat << EOF > $input_file_2

$ntrans
${l1[*]}
${l2[*]}

$ninc
${th_inc[*]}
${ph_inc[*]}

$nvel_raytrace $dvel_raytrace $v0

$star_mode $T_star $R_star $X_star $Y_star $Z_star
$dust_ii_gas
$b_rad $b_tor $b_pol $b_dip
$kmax $extra_info 
$pix_mod $im_pix $im_size
$fits $s_d 
$max_ne_temp $n_depth_rt

EOF

module load intel

cp /priv/c3-astro-data1/PORTAL/bin/main_portal.x ./
./main_portal.x $input_file_1 $output_map $pop $grid $mol $dust 
rm main_portal.x 

#the fits libraries have to be loaded.
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"PATH TO FITSIO LIBRARY"

cp /priv/c3-astro-data1/PORTAL/bin/finish_portal.x ./
./finish_portal.x $input_file_2 $output_map $pop $grid $mol $dust
rm finish_portal.x


