#!/usr/bin/env bash

function psredit_read(){
    datafile=$1
    param_name=$2
    
    param_value=`psredit -c $param_name $datafile 2> /dev/null -q -Q | tr -d ' '`
    echo $param_value 
}

echo "========================================================================="
echo "This script does the following operations."
#echo "1. Convert to PSRFITS"
echo "1. Time collapse"
echo "2. Make bandwidth positive"
#echo "4. Correct source coordinates"
echo "3. Correct the frequency between MJDs 59217 and 59424"
echo "========================================================================="
echo "Ensure that:"
echo "1. It is run on a copy of the original data. It rewrites data files."
echo "2. It is run only once. Otherwise, the frequency will be over-corrected."
echo "========================================================================="
read -n 1 -r -s -p $'Press ENTER to continue, Ctrl+C to abort.\n'

for archive_file in $@
do
    echo
    echo Processing $archive_file ...
    
    # Time collapse, convert to PSRFITS
    echo pam -T -m $archive_file
    pam -T -m $archive_file  
    
    # If bandwidth is negative, reverse channels
    bw=$(psredit_read $archive_file bw)
    if test $bw -le 0
    then
        echo pam --reverse_freqs -m $archive_file
        pam --reverse_freqs -m $archive_file
    fi
    
    # Update the coordinates
    #echo update_coords.sh $archive_file
    #update_coords.sh $archive_file
    
    # Frequency correction
    data_mjd=$(psredit_read $archive_file "int[0]:mjd")
    if test `echo "$data_mjd>=59218 && $data_mjd<59424" | bc -l` == 1
    then
        freq_centre=$(psredit_read $archive_file freq)
        freq_centre_new=`echo $freq_centre + 0.01 | bc -l`
        echo pam -o $freq_centre_new -m $archive_file 2> /dev/null
        pam -o $freq_centre_new -m $archive_file 2> /dev/null
    fi
done

