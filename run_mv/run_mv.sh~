#!/bin/bash
function usage
{
    echo "------------------------- Usage -------------------------"
    echo ""
    echo "      Minimum input is: PPI file and MIPS file with"
    echo "      either a DSD file or a path to DSD. Will run all"
    echo "      necessary scripts to produce a confidence score,"
    echo "      or whatever output is specified. Can optionally"
    echo "      save all intermediate files in output directory,"
    echo "      but will otherwise only save whatever the final"
    echo "      step is that the user specifies."
    echo ""
    echo "-h / --help: Display this help message and exit."
    echo ""
    echo "-i / --inter: Boolean, saves intermediate files."
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
    echo "                    REQUIRED IF annotation or"
    echo "                    permutation files not specified."
    echo ""
    echo "-ad ADJ / --adj ADJ: Adjacency matrix for largest"
    echo "                    connected portion of the network."
    echo "                    Used to start later in the process."
    echo "                    [optional]"
    echo ""
    echo "-t TMAT / --tmat TMAT: Triangular DSD matrix."
    echo ""
    echo "-n NAME / --name NAME: Start of file name that you want"
    echo "                    to be in all output file names."
    echo "                    (e.g. 'biogrid_ppi_')"
    echo "                    REQUIRED"
    echo ""
    echo "-o ODIR / --odir ODIR: Path to output directory."
    echo "                    Defaults to current directory."
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

inter=$false

while :
do
    case $1 in
        -h | --help | -\?)
            usage
            exit 0
            ;;
        -i | --inter)
            inter=$true;;
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

if [[ $mips = "" && $ann = "" && $rand = "" ]]
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
    python2.7 "${dpath}/DSDmain.py" -m 1 -o "${odir}/${name}.DSD1" "${ppi}"
    dsd="${odir}/${name}"
fi

# Find connected PPI
if [[ $cppi = "" ]]
then
    # Run find_connected_list, setting cppi to output
    ./find_connected_list -d "${dsd}" -p "${ppi}" -c "${odir}/${name}_connected.ppi"
    cppi="${odir}/${name}_connected.ppi"
fi

# Find adjacency matrix
if [[ $adj = "" || $oprot = "" ]]
then
    # Run GetAdjacency, setting adj and oprot to output
    adj="Placeholder for adjacency matrix."
    oprot="Placeholdr for ordered protein list."
fi

# Find annotations and create permutation file.
if [[ $ann = "" || $rand = "" ]]
then
    # Run funannotate, setting ann and rand to output
    ann="Placeholder for annotation file."
    rand="Placeholder for permutation file."
fi

# Find the triangular DSD matrix.
if [[ $trimat = "" ]]
then
    # Run make_ordered_trimat, setting trimat to output
    trimat="Placeholder for trimat file."
fi

# Run majority vote, setting mv to output.
mv="Placeholder for MV."

# Run CalculatePerformance
# Output to some file somewhere.
echo "Wrote performance score to w/e file for all 3 levels."