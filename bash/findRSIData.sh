#!/bin/bash
# Purpose: Read Comma Separate CSV File -- Check if RSI data exits and then move cases into a specific folder
# Author: Karoline Kallis
# ------------------------------------------


INPUT="/home/kkallis/RedCap/data/FPRSI_meta_28-Nov-2022.csv"
#INPUT=$1
DataPath="/space/bil-syn01/1/cmig_bil/RSIData/Prostate/UCSD/FPRSIv2/proc/"
DataPath1="/space/bil-syn01/1/cmig_bil/RSIData/Prostate/USCD/FPRSIv2/raw_dicom/"




checkFolder=$1
writeMRN=$2
if [ "$writeMRN" = true ]
then

	[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
	OIFS=$IFS
	IFS=','
	while read PDS_ID FPRSI_ID rest 
	do
       		echo "${FPRSI_ID}" >> "MRN.txt"
	done < $INPUT
	IFS=$OIFS
fi


if [ "$checkFolder" = true ]
then
	cd $DataPath
	
	for dir in */;
	do
		MRN=${dir%*/}
		echo "$MRN" >> "/home/kkallis/Calibration/ProcessedMRN.txt"
		
		if  [[ -n $(find $dir -maxdepth 3 -type f -name "RSI*.mgz") ]]
		then
			if  [ $(grep -w "$MRN" "/home/kkallis/Calibration/MRN.txt") ];
			then
				#echo "$MRN" 

				echo "$MRN" >> "/home/kkallis/Calibration/Processed.txt"
			else
				#echo "$MRN" "Missing Clinical Information!"

				#echo "$MRN" >> "/home/kkallis/Calibration/Processed.txt"
				echo "$MRN" >> "/home/kkallis/Calibration/MissingClinicalInfo.txt"
			fi		
		else
			echo "$MRN" "PROBLEM!!!"
			#echo "$MRN"            "" >> "/home/kkallis/Calibration/RSIDataMissing.txt"
		fi

	done
        #cd $DataPath1

        #for dir in */;
        #do
         #       MRN=${dir%*/}
          #      echo "$MRN"
	#	if ! [ $(grep -w "$MRN" "/home/kkallis/Calibration/Processed.txt") ];
         #       then	
	#		if  [[ -n $(find $dir -maxdepth 2 -type d -name "*RSI*") ]];
	#		then
				#echo $dir
	#			Name=$(find $dir -maxdepth 2 -type d -iname "*RSI*" )
         #                       echo "${MRN} ${Name##*/} " >> "/home/kkallis/Calibration/RSI_DWI.txt"

 	#		else
	#			echo "${MRN}" "No RSI_DWI found!" >>"/home/kkallis/Calibration/RSI_DWI.txt"

				#echo "$MRN" >> "/home/kkallis/Calibration/RSIDataMissing.txt"
	#		fi	
		#echo "$MRN"	"RSI Data Missing"
                #else
                        #echo "$MRN" >>	"/home/kkallis/Calibration/SanityCheck.txt"
                                
         #       fi
	#done
	#find $DataPath -maxdepth 1 -type d | while read dir; do
	#MRN=${dir##*/%*_}
	#echo "${MRN}"


	#if [[ -n $(find $dir -maxdepth 3 -type d -name "RSI*.mgz") ]] 
	#then
	#	echo "${dir##*/}"
	#fi
	#done
fi
