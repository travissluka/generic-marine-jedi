#!/bin/bash
set -eu

#-------------------------------------------------------------------------------
# user configurables that should be set before running
#-------------------------------------------------------------------------------

EXP_DIR=$(readlink -f .) # assume experiment directory is the directory
                         # that this script is in

EXP_START_DATE=20180415Z12  # start date, in YYYYMMDDZHH format
EXP_END_DATE=20180416Z12    # end date,   in YYYYMMDDZHH format

# location of root GenericMarineJEDI source directory
# (not needed by this script, other than the subsequent default IC/landmask locations)
MARINEJEDI_SRC_DIR=$(readlink -f ../)

# location of the GenericMarineJEDI binary files
MARINEJEDI_BIN_DIR=$(readlink -f ../build/bin)

# observation locations, note that %Y %m %d %H will be replaced with the
# date/time of each cycle
OBS_FILE=$EXP_DIR/obs/sst_noaa19_%Y%m%d.nc

# initial background used for the first cycle
IC_FILE=$MARINEJEDI_SRC_DIR/test/data_static/sst_1p0.nc

# land mask and rossby radius file
LANDMASK_FILE=$MARINEJEDI_SRC_DIR/test/data_static/landmask_1p0.nc
ROSSBYRADIUS_FILE=$MARINEJEDI_SRC_DIR/test/data_static/rossby_radius.dat

# temporary working files, that are deleted after each cycle is finished
SCRATCH_DIR=$EXP_DIR/SCRATCH

# MPI command to run will be machine dependent
MPI_PE=8
MPI_CMD="mpirun -n $MPI_PE"

#-------------------------------------------------------------------------------
# END of user configurables
#-------------------------------------------------------------------------------


# you probably dont need to edit anything below here...

# define some other variables used by this script
export OMP_NUM_THREADS=1  # OpenMP isn't working correctly, disable
DA_WINDOW_LEN=24          # assuming a fixed window of 1 day, for now

#================================================================================
#================================================================================
# Start of loop
#================================================================================
#================================================================================
export TIMEFORMAT='%1Rs'
date_fmt="%Y-%m-%dT%H:%M:%SZ"
cycle_avg_count=0
cycle_avg_runtime=0
while true; do
    cycle_start=$(date +%s)

    # determine where the experiment left off and if this is the first cycle
    cycle_status_file=$EXP_DIR/cycle_status
    if [[ ! -e $cycle_status_file ]]; then
        cycle_start_date=$(date -ud "$EXP_START_DATE" +$date_fmt)
        init=1
    else
        cycle_start_date=$(cat $cycle_status_file)
        init=0
    fi

    # are we done with the experiment?
    if [[ $(date -ud "$cycle_start_date" +%Y%m%d%H) -gt \
          $(date -ud "$EXP_END_DATE" +%Y%m%d%H) ]]; then
        echo "Done with the experiment. Last date: $cycle_start_date"
        exit 0
    fi

    # variables that depend on the cycle date
    ANA_DATE=$(date -ud "$cycle_start_date" +$date_fmt)
    DA_WINDOW_HW=$(( DA_WINDOW_LEN/2 ))
    DA_WINDOW_START=$(date -ud "$ANA_DATE - $DA_WINDOW_HW hours" +$date_fmt)
    NEXT_DATE=$(date -ud "$ANA_DATE + $DA_WINDOW_LEN hours" +$date_fmt)
    PREV_DATE_YMDH=$(date -ud "$ANA_DATE - $DA_WINDOW_LEN hours" +"%Y%m%d%H")
    ANA_DATE_YMDH=$(date -ud "$ANA_DATE" +"%Y%m%d%H")
    SCRATCH_DIR_CYCLE=$SCRATCH_DIR/$ANA_DATE_YMDH

    echo ""
    echo "======================================================================"
    echo "analysis date: $ANA_DATE"
    echo "======================================================================"

    # setup working directory
    [[ -d "$SCRATCH_DIR_CYCLE" ]]  && rm -r $SCRATCH_DIR_CYCLE
    mkdir -p $SCRATCH_DIR_CYCLE
    cd $SCRATCH_DIR_CYCLE
    ln -s $LANDMASK_FILE landmask.nc
    ln -s $ROSSBYRADIUS_FILE rossby_radius.dat

    # link in the most recent background state
    if [[ "$init" == 1 ]]; then
        ln -s $IC_FILE bkg.nc
    else
        ln -s $EXP_DIR/ana/ana.$PREV_DATE_YMDH.nc bkg.nc
    fi

    # initialize bump if not already done so
    # TODO, the bump files are dependent on the number of PEs, so if you change PEs, you have
    # to regenerate the bump directory
    BUMP_DIR=$EXP_DIR/bump
    if [[ ! -d "$BUMP_DIR" ]]; then
        echo "Initializing BUMP..."
        mkdir bump
        cp $EXP_DIR/config/{setcorscales,errorcovariance_training}.yaml .        
        
        # generate horizontal length scales, based on rossby radius
        $MPI_CMD $MARINEJEDI_BIN_DIR/genericmarine_setcorscales.x setcorscales.yaml

        # initialize NICAS correlation operator
        $MPI_CMD $MARINEJEDI_BIN_DIR/genericmarine_errorcovariance_training.x errorcovariance_training.yaml

        mv bump $BUMP_DIR
    fi
    ln -s $BUMP_DIR bump

    # link the obserations
    obs_file=$(date -ud "$ANA_DATE" +$OBS_FILE)
    ln -s $obs_file obs.nc    

    # run the var
    mkdir -p obs_out
    cp $EXP_DIR/config/{var,obsop_name_map}.yaml .
    sed -i "s/__DA_WINDOW_START__/${DA_WINDOW_START}/g" var.yaml
    sed -i "s/__ANA_DATE__/${ANA_DATE}/g" var.yaml
    $MPI_CMD $MARINEJEDI_BIN_DIR/genericmarine_var.x var.yaml

    # move the output files
    ana_file=$EXP_DIR/ana/ana.${ANA_DATE_YMDH}.nc
    mkdir -p $(dirname $ana_file)
    mv ana.nc $ana_file

    inc_file=$EXP_DIR/inc/inc.${ANA_DATE_YMDH}.nc
    mkdir -p $(dirname $inc_file)
    mv inc.nc $inc_file

    # moving observation output stats as is
    # TODO: concatenate into a single file per platform
    obs_out_dir=$EXP_DIR/obs_out/${ANA_DATE_YMDH}
    mkdir -p $obs_out_dir
    mv obs_out/*.nc $obs_out_dir

    # done with this day of the cycle, cleanup and prepare for the next cycle
    rm -rf $SCRATCH_DIR_CYCLE
    echo "$NEXT_DATE" > $cycle_status_file

    # update statistics on average cycle runtime
    cycle_end=$(date +%s)
    cycle_runtime=$((cycle_end-cycle_start ))
    cycle_avg_runtime=$(((cycle_avg_runtime*cycle_avg_count+cycle_runtime)/(cycle_avg_count+1) ))
    cycle_avg_count=$((cycle_avg_count+1))
    echo "Cycle runtime:     $cycle_runtime seconds"
    echo "Cycle runtime avg: $cycle_avg_runtime seconds"
done