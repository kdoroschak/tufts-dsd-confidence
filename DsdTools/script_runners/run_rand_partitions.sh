#!/bin/bash -e
function usage {
	echo "                        USAGE"
	echo "---------------------------------------------------------------------"
	echo ""
	echo "-h / --help                       Display this message and quit."
	echo ""
	echo "-ap ANNPATH / --annpath ANNPATH   Pathway to folder containing"
	echo "                                  annotations, including file names"
	echo "                                  up to but excluding '.mips..."
	echo ""
	echo "-mvp MVP / --mvpath MVP           Pathway to majorityvote directory."
	echo ""
	echo "-o ODIR / --odir ODIR             Output directory."
	echo ""
	echo "-n NAME / --name NAME             Name for output files."
	echo ""
	echo "-p PART / --part PART             Pathway to folder containing"
	echo "                                  subdirectories of mips 1-3"
	echo "                                  random partition files."
	echo ""
	echo "-a ADJ / --adj ADJ                Adjacency matrix."
	echo ""
	echo "-d DSD / --dsd DSD                DSD triangular matrix."
	echo "---------------------------------------------------------------------"

}

function quit {
	exit
}

trap quit SIGINT

while :
do
    case $1 in
        -h | --help | -\?)
            usage
            exit 0
            ;;
		-ap | --annpath)
			if [ $2 ]
			then
				annpath=$2
				shift 2
			fi
			;;
		-mvp | --mvpath)
			if [ $2 ]
			then
				mvp=$2
				shift 2
			fi
			;;
		-o | --odir)
			if [ $2 ]
			then
				odir=$2
				shift 2
			fi
			;;
		-n | --name)
			if [ $2 ]
			then
				name=$2
				shift 2
			fi
			;;
		-p | --part)
			if [ $2 ]
			then
				part=$2
				shift 2
			fi
			;;
		-a | --adj)
			if [ $2 ]
			then
				adj=$2
				shift 2
			fi
			;;
		-d | --dsd)
			if [ $2 ]
			then
				dsd=$2
				shift 2
			fi
			;;
        --)
            shift
            break
            ;;
        -*)
            echo "Warning: Unknown option: $1. Ignoring."
            echo "-h / --help for usage."
            shift
            ;;
        *)
            break
            ;;
    esac
done

failed=false

if [[ $mvp = "" ]]
then
	echo "Error: No majority vote path specified."
	failed=true
fi

if [[ $annpath == "" ]]
then
	echo "Error: No pathway to annotations specified."
	failed=true
fi

if [[ $name = "" ]]
then
    echo "Error: No name specified."
    failed=true
fi

if [[ $odir = "" ]]
then
	echo "Error: No output directory specified."
	failed=true
fi

if [[ $part = "" ]]
then
	echo "Error: No partition directory specified."
	failed=true
fi

if [[ $adj = "" ]]
then
	echo "Error: No adjacency matrix specified."
	failed=true
fi

if [[ $dsd = "" ]]
then
    echo "Error: No DSD triangular matrix specified."
    failed=true
fi

if $failed
then
    echo "FAILED"
    echo "-h/--help for usage."
    exit
fi

ann1="${annpath}.mips1.ann"
ann2="${annpath}.mips2.ann"
ann3="${annpath}.mips3.ann"

for i in $(seq -f "%02g" 0 99)
do
	part1="${part}/rand_partition_test_mips1/rand_partition_test-$i.mips1.ann.rand"
	part2="${part}/rand_partition_test_mips2/rand_partition_test-$i.mips2.ann.rand"
	part3="${part}/rand_partition_test_mips3/rand_partition_test-$i.mips3.ann.rand"
	subdir="random_partition_test_$i"
	out="${odir}/${subdir}"
	if [ ! -d "${odir}/${subdir}" ]
	then
		mkdir "${odir}/${subdir}"
	fi

	echo " ------- Running majorityvote --------"

	python2.7 "${mvp}/DSDmv.py" -l $ann1 -r $part1 -k 2 -o "_${name}-$i.mips1.mv" -t 10 -m 2 -d $dsd $adj
	python2.7 "${mvp}/DSDmv.py" -l $ann2 -r $part2 -k 2 -o "_${name}-$i.mips2.mv" -t 10 -m 2 -d $dsd $adj
	python2.7 "${mvp}/DSDmv.py" -l $ann3 -r $part3 -k 2 -o "_${name}-$i.mips3.mv" -t 10 -m 2 -d $dsd $adj

	mv DSDWeighted_${name}-$i.mips* $out

	echo " ------ Calculating performance ------"

	python2.7 "${mvp}/CalculatePerformance.py" $ann1 "${out}/DSDWeighted_${name}-$i.mips1.mv" > "${out}.perf"
	python2.7 "${mvp}/CalculatePerformance.py" $ann2 "${out}/DSDWeighted_${name}-$i.mips2.mv" > "${out}.perf"
	python2.7 "${mvp}/CalculatePerformance.py" $ann3 "${out}/DSDWeighted_${name}-$i.mips3.mv" > "${out}.perf"

	chmod 664 ${odir}/${subdir}/*
	chgrp bcb ${odir}/${subdir}/*
	chmod 775 ${odir}/${subdir}
	chgrp bcb ${odir}/${subdir}

	echo " ---- Performacne scores written -----"

done
