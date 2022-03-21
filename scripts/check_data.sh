# A script to check that the data in GLOBAL_MULTIYEAR and GLOBAL_REANALYSIS are
# the same.

# Only a subset is considered - the daily averages.
# One year at a time.
# Usage: bash check_data.sh YYYY

# Nancy Soontiens Feb 2022

multiyear=/data/cmems/my.cmems-du.eu/GLOBAL_MULTIYEAR_PHY_001_030/cmems_mod_glo_phy_my_0.083_P1D-m/
reanalysis=/data/cmems/my.cmems-du.eu/GLOBAL_REANALYSIS_PHY_001_030/global-reanalysis-phy-001-030-daily

year=$1
sd=${year}0101
start=$(date -d $sd +%Y%m%d)
end=$(date -d "$start + 1 year" +%Y%m%d)

logfile=${year}_check.out
echo "Comparing MULTIYEAR and REANALYSIS for each day in ${year} using: " > ${logfile}
echo "diff multiyear_file reanalysis_file" >> ${logfile}
echo "A blank line means no difference" >> ${logfile}


while [[ $start != $end ]]; do
    
    month=$(date -d"$start" +"%m")

    multiyear_file=${multiyear}/${year}/${month}/*${start}_*.nc
    reanalysis_file=${reanalysis}/${year}/${month}/*${start}_*.nc
    echo "" >> ${logfile}
    echo "${start}" >> ${logfile}
    if [ ! -f ${multiyear_file} ]; then
	echo "MULTIYEAR file for ${start} does not exist" >> ${logfile}
	echo ${multiyear_file} >> multiyear_missingfiles
    fi
    if [ ! -f ${reanalysis_file} ]; then
	echo "REANALYSIS file for ${start} does not exist" >> ${logfile}
        echo ${reanalysis_file} >> reanalysis_missingfiles
    fi
    if [ -f ${reanalysis_file} ] && [ -f ${multiyear_file}  ]; then
	echo "diff check:" >> ${logfile}
	diff ${multiyear_file} ${reanalysis_file} >> ${logfile}
	if  ! ncdump -h ${multiyear_file} ; then
	    echo ${multiyear_file} >> multiyear_badfiles
	fi
	if  ! ncdump -h ${reanalysis_file} ; then
            echo ${reanalysis_file} >> reanalysis_badfiles
        fi

    fi
    start=$(date -d"$start + 1 day" +"%Y%m%d")
done
    

