#!/bin/bash 

#PBS -N run_l96_exps
#PBS -e run_l96_exps.err
#PBS -o run_l96_exps.out
#PBS -l select=XX:ncpus=XX:mpiprocs=XX
#PBS -l walltime=XX:XX:XX
#PBS -A XXXXXXXXX
#PBS -q regular
#PBS -m abe
#PBS -M gharamti@ucar.edu


##############################################################
# Prior to executing this script, one needs to compile DART: # 
# ./quickbuild.sh                                            #
# The script below does a lot more experiments than what's   #
# shown in the paper. Use with caution!                      #
##############################################################


export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

WORKDIR=/glade/work/$USER/DART/dart_rd/models/lorenz_96/work

cd $WORKDIR

nproc=10
filter=1
expall=true

if $expall ; then
   suff=all
   ./create_obs_sequence      < ../obs_input/identity_obs.input
   ./create_fixed_network_seq < ../obs_input/fixednet.input
   ./perfect_model_obs
else
   suff=half
   ./create_obs_sequence      < ../obs_input/every_other_obs.input
   ./create_fixed_network_seq < ../obs_input/fixednet.input
   ./perfect_model_obs
fi

if [[ $filter == 1 ]]; then 
   filter_name=EAKF_$suff
   increment=.false.
elif [[ $filter == 2 ]]; then 
   filter_name=EnKF_$suff
   increment=.true.
fi

DIAGDIR=/glade/scratch/$USER/rd_enkf/L96/${filter_name}

mkdir -p $DIAGDIR

# Copy true state to diagnostic directory
cp ${WORKDIR}/true_state.nc ${DIAGDIR}

echo -e " "
echo -e "Selected Filter: ${filter_name}; 'sort_obs_inc' is set to ${increment}"
echo -e "Script is starting ..."
echo -e " "

# Select the right filter kind: 
sed -i "s/.*filter_kind.*/   filter_kind                     = ${filter},/" input.nml
sed -i "s/.*sort_obs_inc.*/   sort_obs_inc                    = ${increment},/" input.nml

do_enss=true
do_perf=false # 20 members
do_impr=false # 20 members
do_rtps=true
do_obsn=false
do_obse=false
do_obst=false

# 1. Ensemble Size: 
# *****************
if $do_enss ; then 

   edir=$DIAGDIR/enssens_up
   mkdir -p $edir

   # No Model Error
   sed -i "s/.*forcing.*/   forcing           = 8,/" $WORKDIR/input.nml

   # localization
   sed -i "s/.*cutoff .*/   cutoff                          = 0.10,/" $WORKDIR/input.nml

   # Randomization every obs
   sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .true.,/" $WORKDIR/input.nml

   # Inflation flavor
   sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       0,/" $WORKDIR/input.nml

   ens=($(seq -w 5 5 80))
   Ne=`echo "(${#ens[@]} - 1)" | bc -l`

   alpha=( 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 )
   Na=`echo "(${#alpha[@]} - 1)" | bc -l`

   for jj in `seq 0 $Na`; do  
       sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" $WORKDIR/input.nml
   
       for ii in `seq 0 $Ne`; do        
   
           sed -i "s/.*ens_size .*/   ens_size                     = ${ens[$ii]},/" $WORKDIR/input.nml

           # i- No inflation:
           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" $WORKDIR/input.nml

           mpirun -np $nproc ./filter ; wait  
           mv preassim.nc ${edir}/preassim_alp${alpha[$jj]}_ens${ens[$ii]}.nc
           mv analysis.nc ${edir}/analysis_alp${alpha[$jj]}_ens${ens[$ii]}.nc

           # ii- Add inflation: 
           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.6,                     0.0,/" $WORKDIR/input.nml

           mpirun -np $nproc ./filter ; wait
           mv preassim.nc ${edir}/preassim_alp${alpha[$jj]}_ens${ens[$ii]}_inf.nc
           mv analysis.nc ${edir}/analysis_alp${alpha[$jj]}_ens${ens[$ii]}_inf.nc

       done
   done
fi

# 2. No Model Errors: 
# *******************

if $do_perf ; then 

   pdir=$DIAGDIR/perfect
   mkdir -p $pdir

   # No model errors
   sed -i "s/.*forcing.*/   forcing           = 8,/" input.nml

   # Randomization every obs
   sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .true.,/" input.nml

   # Inflation flavor
   sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       0,/" input.nml  

   # Fix the ensemble size
   sed -i "s/.*ens_size .*/   ens_size                     = 10,/" input.nml
 
   cut=($(seq -w 0.10 0.02 0.30)) #( 0.1 0.2 0.3 0.4 ) 
   Nl=`echo "(${#cut[@]} - 1)" | bc -l`
   
   alpha=( 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 )
   Na=`echo "(${#alpha[@]} - 1)" | bc -l`
   
   for jj in `seq 0 $Na`; do 
       sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" input.nml
   
       for ii in `seq 0 $Nl`; do    
   
           sed -i "s/.*cutoff .*/   cutoff                          = ${cut[$ii]},/" input.nml    
  
           # i- No inflation:
           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" input.nml           
 
           mpirun -np $nproc ./filter ; wait    
           mv preassim.nc ${pdir}/preassim_alp${alpha[$jj]}_loc${cut[$ii]}.nc    
           mv analysis.nc ${pdir}/analysis_alp${alpha[$jj]}_loc${cut[$ii]}.nc    

           # ii- Add inflation: 
           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.6,                     0.0,/" input.nml

           mpirun -np $nproc ./filter ; wait    
           mv preassim.nc ${pdir}/preassim_alp${alpha[$jj]}_loc${cut[$ii]}_inf.nc    
           mv analysis.nc ${pdir}/analysis_alp${alpha[$jj]}_loc${cut[$ii]}_inf.nc  

       done
   done
fi

# 3. With Model Errors: 
# *********************

if $do_impr ; then 

   idir=$DIAGDIR/imperfect
   mkdir -p $idir

   # Fix the localization radius
   sed -i "s/.*cutoff .*/   cutoff                          = 0.10,/" input.nml

   # Fix the ensemble size
   sed -i "s/.*ens_size .*/   ens_size                     = 20,/" input.nml

   # Inflation flavor
   sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       0,/" input.nml

   alpha=( 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 )
   #alpha=( 0.4 0.5 0.6 0.7 0.8 )
   Na=`echo "(${#alpha[@]} - 1)" | bc -l`

   forcing=( 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 )
   Nf=`echo "(${#forcing[@]} - 1)" | bc -l`

   # (1) Randomize every obs
   sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .true.,/" input.nml

   for jj in `seq 0 $Na`; do
       sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" input.nml

       for ii in `seq 0 $Nf`; do 
           sed -i "s/.*forcing.*/   forcing           = ${forcing[$ii]},/" input.nml

           # i- No inflation
           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" input.nml

           mpirun -np $nproc ./filter ; wait 
           mv preassim.nc ${idir}/preassim_alp${alpha[$jj]}_for${forcing[$ii]}_T.nc
           mv analysis.nc ${idir}/analysis_alp${alpha[$jj]}_for${forcing[$ii]}_T.nc

           # ii- Add inflation
           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.6,                     0.0,/" input.nml

           mpirun -np $nproc ./filter ; wait 
           mv preassim.nc ${idir}/preassim_alp${alpha[$jj]}_for${forcing[$ii]}_inf_T.nc
           mv analysis.nc ${idir}/analysis_alp${alpha[$jj]}_for${forcing[$ii]}_inf_T.nc
       done
   done

#   # (2) Randomize only once
#   sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .false.,/" input.nml
#
#   for jj in `seq 0 $Na`; do
#       sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" input.nml
#
#       for ii in `seq 0 $Nf`; do 
#           sed -i "s/.*forcing.*/   forcing           = ${forcing[$ii]},/" input.nml
#
#           # i- No inflation
#           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" input.nml
#
#           mpirun -np $nproc ./filter ; wait 
#           mv preassim.nc ${idir}/preassim_alp${alpha[$jj]}_for${forcing[$ii]}_F.nc
#           mv analysis.nc ${idir}/analysis_alp${alpha[$jj]}_for${forcing[$ii]}_F.nc
#
#           # ii- Add inflation
#           sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.6,                     0.0,/" input.nml
#
#           mpirun -np $nproc ./filter ; wait 
#           mv preassim.nc ${idir}/preassim_alp${alpha[$jj]}_for${forcing[$ii]}_inf_F.nc
#           mv analysis.nc ${idir}/analysis_alp${alpha[$jj]}_for${forcing[$ii]}_inf_F.nc
#       done
#   done

fi


# 4. Compare to RTPS: 
# *******************

if $do_rtps ; then

   rdir=$DIAGDIR/rtps_N20_F3_up
   mkdir -p $rdir

   # Fix the localization radius
   sed -i "s/.*cutoff .*/   cutoff                          = 0.10,/" input.nml

   # Fix the ensemble size
   sed -i "s/.*ens_size .*/   ens_size                     = 20,/" input.nml

   # inflation
   sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       4,/" input.nml
   sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" input.nml

   # No RD
   alpha=0.0
   sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha},/" input.nml
   sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .true.,/" input.nml   

   # ensemble size
   ens=($(seq -w 10 10 80))
   Ne=`echo "(${#ens[@]} - 1)" | bc -l`

   # RTPS weighting factor
   weight=($(seq -w 0.00 0.10 1.01))
   Nw=`echo "(${#weight[@]} - 1)" | bc -l`

   # model error
   sed -i "s/.*forcing.*/   forcing           = 3,/" $WORKDIR/input.nml

   # RTPS
   for jj in `seq 0 $Nw`; do
       sed -i "s/.*inf_initial .*/   inf_initial                 = 1.0,                     ${weight[$jj]},/" input.nml

       for ii in `seq 0 $Ne`; do
           sed -i "s/.*ens_size .*/   ens_size                     = ${ens[$ii]},/" input.nml

           mpirun -np $nproc ./filter ; wait
           mv preassim.nc ${rdir}/preassim_ens${ens[$ii]}_rtps${weight[$jj]}.nc
           mv analysis.nc ${rdir}/analysis_ens${ens[$ii]}_rtps${weight[$jj]}.nc
      done
   done

   alpha=( 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 )
   Na=`echo "(${#alpha[@]} - 1)" | bc -l`

   # no inflation
   sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       0,/" input.nml

   # RD
   for jj in `seq 0 $Na`; do
       sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" input.nml

       for ii in `seq 0 $Ne`; do
           sed -i "s/.*ens_size .*/   ens_size                     = ${ens[$ii]},/" input.nml

           mpirun -np $nproc ./filter ; wait
           mv preassim.nc ${rdir}/preassim_ens${ens[$ii]}_alp${alpha[$jj]}.nc
           mv analysis.nc ${rdir}/analysis_ens${ens[$ii]}_alp${alpha[$jj]}.nc
      done
   done
fi
 

# 5. Obs Network: 
# ***************

if $do_obsn ; then 

  ndir=$DIAGDIR/obsnetwork_F8_N10
  mkdir -p $ndir

  # fix localization 
  sed -i "s/.*cutoff .*/   cutoff                          = 0.10,/" input.nml

  # fix forcing: no model errors
  sed -i "s/.*forcing.*/   forcing           = 8,/" input.nml

  # fix alpha 
  alpha=( 0.00 0.05 0.10 0.15 0.20 0.25 0.30 )
  Na=`echo "(${#alpha[@]} - 1)" | bc -l`

  # fix the ensemble size
  sed -i "s/.*ens_size .*/   ens_size                     = 10,/" input.nml

  # randomization
  sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .true.,/" input.nml

  # trun off RTPS
  sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       0,/" input.nml

  # 6 networks 
  networks=( 'every_tenth' 'mid_five' 'every_fifth' 'last_ten' 'every_fourth' 'every_third' 'every_other' 'identity' )
  Nl=`echo "(${#networks[@]} - 1)" | bc -l`  

  for ii in `seq 0 $Nl`; do
      net=${networks[$ii]}
    
      # prepare obs, run pmo
      ./create_obs_sequence < ../obs_input/${net}_obs.input
      ./create_fixed_network_seq < ../obs_input/fixednet.input
      ./perfect_model_obs ; wait

      for jj in `seq 0 $Na`; do
          # run filter
          sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" input.nml
 
          # i- No inflation
          sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" input.nml
    
          mpirun -np $nproc ./filter ; wait 
          mv preassim.nc ${ndir}/preassim_alp${alpha[$jj]}_obs${ii}.nc
          mv analysis.nc ${ndir}/analysis_alp${alpha[$jj]}_obs${ii}.nc
    
          # ii- Add inflation
          sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.6,                     0.0,/" input.nml
    
          mpirun -np $nproc ./filter ; wait 
          mv preassim.nc ${ndir}/preassim_alp${alpha[$jj]}_obs${ii}_inf.nc
          mv analysis.nc ${ndir}/analysis_alp${alpha[$jj]}_obs${ii}_inf.nc
      done 
  done    

fi


# 6. Obs Error Variance: 
# **********************

if $do_obse ; then

  ndir=$DIAGDIR/obserror_N10
  mkdir -p $ndir

  # fix localization 
  sed -i "s/.*cutoff .*/   cutoff                          = 0.10,/" input.nml

  # fix forcing: no model errors
  sed -i "s/.*forcing.*/   forcing           = 8,/" input.nml

  # fix alpha 
  alpha=( 0.00 0.05 0.10 0.15 0.20 0.25 0.30 )
  Na=`echo "(${#alpha[@]} - 1)" | bc -l`

  # fix the ensemble size
  sed -i "s/.*ens_size .*/   ens_size                     = 10,/" input.nml

  # Only prior inflation 
  sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       0,/" input.nml

  # randomization
  sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .true.,/" input.nml

  # 6 networks 
  networks=( 'identity_1' 'identity' 'identity_4' 'identity_6' 'identity_8' 'identity_10' )
  Nl=`echo "(${#networks[@]} - 1)" | bc -l`

  for ii in `seq 0 $Nl`; do
      net=${networks[$ii]}

      # prepare obs, run pmo
      ./create_obs_sequence < ../obs_input/${net}_obs.input
      ./create_fixed_network_seq < ../obs_input/fixednet.input
      ./perfect_model_obs ; wait

      for jj in `seq 0 $Na`; do
          # run filter
          sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" input.nml

          # i- No inflation
          sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" input.nml

          mpirun -np $nproc ./filter ; wait
          mv preassim.nc ${ndir}/preassim_alp${alpha[$jj]}_obs${ii}.nc
          mv analysis.nc ${ndir}/analysis_alp${alpha[$jj]}_obs${ii}.nc

          # ii- Add inflation
          sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.6,                     0.0,/" input.nml

          mpirun -np $nproc ./filter ; wait
          mv preassim.nc ${ndir}/preassim_alp${alpha[$jj]}_obs${ii}_inf.nc
          mv analysis.nc ${ndir}/analysis_alp${alpha[$jj]}_obs${ii}_inf.nc
      done
  done

fi


# 7. Obs Frequency: 
# *****************

if $do_obst ; then

  ndir=$DIAGDIR/obfreq_N20_F10
  mkdir -p $ndir

  # fix localization 
  sed -i "s/.*cutoff .*/   cutoff                          = 0.10,/" input.nml

  # fix forcing: no model errors
  sed -i "s/.*forcing.*/   forcing           = 8,/" input.nml

  # fix alpha 
  alpha=( 0.00 0.05 0.10 0.15 0.20 0.25 0.30 )
  Na=`echo "(${#alpha[@]} - 1)" | bc -l`

  # fix the ensemble size
  sed -i "s/.*ens_size .*/   ens_size                     = 20,/" input.nml

  # Only prior inflation 
  sed -i "s/.*inf_flavor .*/   inf_flavor                  = 5,                       0,/" input.nml

  # randomization
  sed -i "s/.*rd_every_ob.*/   rd_every_ob                     = .true.,/" input.nml

  # 6 networks 
  networks=( '1hr' '2hr' '3hr' '4hr' '5hr' '6hr' '12hr' '18hr' )
  Nl=`echo "(${#networks[@]} - 1)" | bc -l`

  for ii in `seq 0 $Nl`; do
      net=${networks[$ii]}

      # prepare obs, run pmo
      ./create_obs_sequence < ../obs_input/identity_obs.input
      ./create_fixed_network_seq < ../obs_input/fixednet_${net}.input
      ./perfect_model_obs ; wait

      for jj in `seq 0 $Na`; do
          # run filter
          sed -i "s/.*rd_factor.*/   rd_factor                       = ${alpha[$jj]},/" input.nml

          # i- No inflation
          sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.0,                     0.0,/" input.nml

          mpirun -np $nproc ./filter ; wait
          mv preassim.nc ${ndir}/preassim_alp${alpha[$jj]}_obs${ii}.nc
          mv analysis.nc ${ndir}/analysis_alp${alpha[$jj]}_obs${ii}.nc

          # ii- Add inflation
          sed -i "s/.*inf_sd_initial .*/   inf_sd_initial              = 0.6,                     0.0,/" input.nml

          mpirun -np $nproc ./filter ; wait
          mv preassim.nc ${ndir}/preassim_alp${alpha[$jj]}_obs${ii}_inf.nc
          mv analysis.nc ${ndir}/analysis_alp${alpha[$jj]}_obs${ii}_inf.nc
      done
  done

fi
