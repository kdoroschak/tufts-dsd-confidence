#!/bin/bash -e
function usage
{
    echo "------------------------- Usage -------------------------"
    echo ""
    echo "      Minimum input is: PPI file and MIPS file with"
    echo "      either a DSD file or a path to DSD. Will run all"
    echo "      necessary scripts to produce a confidence score."
    echo ""
    echo "-h / --help: Display this help message and exit."
    echo ""
    echo "-p PPI / --ppi PPI: PPI is a PPI input file. This should"
    echo "                    be the file that the DSD was run on."
    echo "                    REQUIRED"
    echo ""
    echo "-d DSD / --dsd DSD: DSD is a DSD1 input file."
    echo "                    If unspecified, will run DSD on PPI."
    echo "                    REQUIRED IF dpath not specified."
    echo ""
    echo "-dp DPATH / --dpath DPATH: Path to DSD directory."
    echo "                    REQUIRED IF dsd not specified."
    echo ""
    echo "-mvp MVPATH / --mvpath MVPATH: Path to majorityvote."
    echo "                    REQUIRED."
    echo ""
    echo "-m MIPS / --mips MIPS: MIPS annotations input file."
    echo "                    REQUIRED IF no path for annotations"
    echo "                    and random permuations directory."
    echo ""
    echo "-ap ANNPATH / --annpath ANNPATH: Path to annotations and"
    echo "                    random permutations directory. Also"
    echo "                    includes everything in the file"
    echo "                    names before '.mipsX.ann[.rand]'."
    echo "                    File names must be identical before"
    echo "                    '.mipsX.ann[.rand]' and must have"
    echo "                    suffixes in that format."
    echo "                    REQUIRED IF no MIPS file specified."
    echo ""
    echo "-ad ADJ / --adj ADJ: Adjacency matrix for largest"
    echo "                    connected portion of the network."
    echo "                    Used to start later in the process."
    echo "                    [optional]"
    echo ""
    echo "-t TMAT / --tmat TMAT: Triangular DSD matrix."
    echo "                    [optional]"
    echo ""
    echo "-n NAME / --name NAME: Start of file name that you want"
    echo "                    to be in all output file names."
    echo "                    (e.g. 'biogrid_ppi_')"
    echo "                    REQUIRED"
    echo ""
    echo "-o ODIR / --odir ODIR: Path to output directory."
    echo "                    Defaults to current directory."
    echo "                    [optional]"
    echo ""
    echo "-cp CPPI / --cppi CPPI: Connected protein list."
    echo "                    Used to start later in the process."
    echo "                    [optional]"
    echo ""
    echo "-op OPROT / --oprot OPROT: Ordered protein list."
    echo "                    Used to start later in the process."
    echo "                    [optional]"
    echo ""
    echo "---------------------------------------------------------"
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
        -p | --ppi)
            if [ $2 ]
            then
                ppi=$2
                shift 2
            else
                echo "No argument for -p/--ppi. FATAL."
                echo "-h/--help for usage."
                exit 0
            fi
            ;;
        -d | --dsd)
            if [ $2 ]
            then
                dsd=$2
                shift 2
            else
                echo "No argument for -d/--dsd. Ignoring."
                shift
            fi
            ;;
        -ap | --annpath)
            if [ $2 ]
            then
                annpath=$2
                shift 2
            else
                echo "No argument for -ap/--annpath. Ignoring."
                shift
            fi
            ;;
        -dp | --dpath)
            if [ $2 ]
            then
                dpath=$2
                shift 2
            else
                echo "No argument for -dp/-dpath. Ignoring."
                shift
            fi
            ;;
        -mvp | --mvpath)
            if [ $2 ]
            then
                mvpath=$2
                shift 2
            else
                echo "No argument for -mvp/--mvpath. FATAL."
                echo "-h/--help for usage."
                exit 0
            fi
            ;;
        -m | --mips)
            if [ $2 ]
            then
                mips=$2
                shift 2
            else
                echo "No argument for -m/--mips. Ignoring."
                shift
            fi
            ;;
        -ad | --adj)
            if [ $2 ]
            then
                adj=$2
                shift 2
            else
                echo "No argument for -ad/--adj. Ignoring."
                shift
            fi
            ;;
        -t | --tmat)
            if [ $2 ]
            then
                trimat=$2
                shift 2
            else
                echo "No argument for -t/--tmat. Ignoring."
                shift
            fi
            ;;
        -n | --name)
            if [ $2 ]
            then
                name=$2
                shift 2
            else
                echo "No argument for -n/--name. FATAL."
                echo "-h/--help for usage."
                exit 0
            fi
            ;;
        -o | --odir)
            if [ $2 ]
            then
                odir=$2
                shift 2
            else
                echo "No argument for -o/--odir. Ignoring."
                shift
            fi
            ;;
        -cp | --cppi)
            if [ $2 ]
            then
                cppi=$2
                shift 2
            else
                echo "No argument for -cp/--cppi. Ignoring."
                shift
            fi
            ;;
        -op | --oprot)
            if [ $2 ]
            then
                oprot=$2
                shift 2
            else
                echo "No argument for -op/--oprot. Ignoring."
                shift
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

if [[ $ppi = "" ]]
then
    echo "Error: No PPI file specified."
    failed=true
fi

if [[ $name = "" ]]
then
    echo "Error: No name specified."
    failed=true
fi

if [[ $mips = "" && $annpath = "" ]]
then
    echo "Error: MIPS not specified."
    failed=true
fi 

if [[ $dsd = "" && $dpath = "" ]]
then
    echo "Error: Neither DSD nor DPATH specified."
    failed=true
fi

if $failed
then
    echo "FAILED"
    echo "-h/--help for usage."
    exit 0
fi

if [[ $odir = "" ]]
then
    odir="."
fi


# Path:
#   DSD
#   Find connected list
#   Get adjacency
#   Annotations/Randomized
#   Trimat
#   MV
#   CalculatePerformance


# Find DSD
if [[ $dsd = "" ]]
then
    # Run DSD, setting $dsd to the output.
    dsd="${odir}/${name}.DSD1"
    python2.7 ${dpath}/DSDmain.py -m 1 -o ${dsd} ${ppi}
fi

# Find connected PPI
if [[ $cppi = "" ]]
then
    # Run find_connected_list, setting cppi to output
    cppi="${odir}/${name}_connected.ppi"
    ./find_connected_list -d ${dsd} -p ${ppi} -c ${cppi}
fi

# Find adjacency matrix
if [[ $adj = "" || $oprot = "" ]]
then
    # Run GetAdjacency, setting adj and oprot to output
    adj="${odir}/${name}_connected.adj"
    oprot="${odir}/${name}_connected_ordered.protein"
    python2.7 getadjacency.py -p ${cppi} -m ${adj} -o ${oprot}
fi

# Run funannotate, setting ann and rand to output
if [[ $annpath = "" ]]
then
    annpath="${odir}/${name}_connected_ordered"
    ./funannotate.py ${mips} ${oprot} -l 1 -o ${annpath}.mips1.ann -p ${annpath}.mips1.ann.rand
    ./funannotate.py ${mips} ${oprot} -l 2 -o ${annpath}.mips2.ann -p ${annpath}.mips2.ann.rand
    ./funannotate.py ${mips} ${oprot} -l 3 -o ${annpath}.mips3.ann -p ${annpath}.mips3.ann.rand
fi

# Find the triangular DSD matrix.
if [[ $trimat = "" ]]
then
    # Run make_ordered_trimat, setting trimat to output
    trimat="${odir}/${name}.dsd1.trimat"
    ./make_ordered_trimat -d ${dsd} -p ${oprot} -m ${trimat}
fi

# Run majority vote, setting mv to output.
python2.7 ${mvpath}/DSDmv.py -l ${annpath}.mips1.ann -r ${annpath}.mips1.ann.rand -k 2 -o _${name}.mips1.mv -d ${trimat} -t 10 -m 2 ${adj}
python2.7 ${mvpath}/DSDmv.py -l ${annpath}.mips2.ann -r ${annpath}.mips2.ann.rand -k 2 -o _${name}.mips2.mv -d ${trimat} -t 10 -m 2 ${adj}
python2.7 ${mvpath}/DSDmv.py -l ${annpath}.mips3.ann -r ${annpath}.mips3.ann.rand -k 2 -o _${name}.mips3.mv -d ${trimat} -t 10 -m 2 ${adj}
mv DSDWeighted_${name}.* ${odir}

# Run CalculatePerformance
# Output to some file somewhere.
python2.7 ${mvpath}/CalculatePerformance.py ${annpath}.mips1.ann ${odir}/DSDWeighted_${name}.mips1.mv > "${odir}/${name}.mips1.perf"
python2.7 ${mvpath}/CalculatePerformance.py ${annpath}.mips2.ann ${odir}/DSDWeighted_${name}.mips2.mv > "${odir}/${name}.mips2.perf"
python2.7 ${mvpath}/CalculatePerformance.py ${annpath}.mips3.ann ${odir}/DSDWeighted_${name}.mips3.mv > "${odir}/${name}.mips3.perf"

echo "Wrote performance score to output directory for all 3 MIPS levels."
