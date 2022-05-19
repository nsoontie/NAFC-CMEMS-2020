# Script to make monthly means of Heather's subsetted data
src_dir=/ocean/handres/glorys12
out_dir=/home/soontiensn/data/cmems2020-oceanstatereport/data/netcdf/GLORYSv12/monthly-means

year=$1
    for month in {01..12}; do
	infile=${src_dir}/Allvars_glorys12_${year}_daily${month}.nc
	echo $infile
	if [[ -f $infile ]]; then
	    outfile=${out_dir}/AllVars_glorys12_${year}_mean${month}.nc
            cdo monmean $infile $outfile
	fi
    done
cdo -seasmean -select,season=JJA -mergetime ${out_dir}/AllVars*${year}* ${out_dir}/${year}_JJA.nc

