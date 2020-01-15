###############################################################################
###############################################################################
#############################INPUT PARAMETERS##################################
###############################################################################
###############################################################################

#input file names

#input from LIME for simulation
pop='/c3se/users/lankhaar/Vera/pol_files/lime/agb/agb_12CO-v2_large.pop'
grid='/c3se/users/lankhaar/Vera/pol_files/lime/agb/agb_12CO-v2_large.vtk'
dust='~/pol_files/jena_thin_e6.tab'
mol='~/pol_files/CO_coro_v2.dat'

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

#integration parameters
nth=40
nph=40
nvel=10

#parameter-length parameters
n_depth=10000
n_depth_rt=10000
max_ne_temp=10000

#maximum irreducible tensor element (>6 is advised)
kmax=6

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
im_pix=100	#ray-tracing resolution 
im_size='0.7'   #image size as a fraction of the simulation size
		#cannot be higher than 1/sqrt(2)

#pix_mod=0 ray-traces the solution at all of the node-points projected on the
#plane of the sky. The output files then will be accompanied 

delta_pix="2.7778e-4"

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
$b_rad $b_tor $b_pol $b_dip
$kmax $extra_info 
$max_ne_temp $n_depth

EOF

cat << EOF > $input_file_2

$ntrans
${l1[*]}
${l2[*]}

$ninc
${th_inc[*]}
${ph_inc[*]}

$nvel_raytrace $dvel_raytrace $v0

$star_mode $T_star $R_star $X_star $Y_star $Z_star
$b_rad $b_tor $b_pol $b_dip
$kmax $extra_info 
$pix_mod $im_pix $im_size $s_d
$max_ne_temp $n_depth_rt

EOF

module load intel

cp /c3se/users/lankhaar/Hebbe/PORTAL/bin/main_portal.x ./
./main_portal.x $input_file_1 $output_map $pop $grid $mol $dust 
rm main_portal.x 

cp /c3se/users/lankhaar/Hebbe/PORTAL/bin/finish_portal.x ./
./finish_portal.x $input_file_2 $output_map $pop $grid $mol $dust
rm finish_portal.x


