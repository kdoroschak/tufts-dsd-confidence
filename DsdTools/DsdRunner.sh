#!/bin/bash

# biogrid_DSD_runner.sh
# Run a batch of files through DSD

# See below the script for instructions.

for i in $(seq -f "%03g" START STOP)
do
	output_file=/path/to/output/DSD/file/example/filename-$i
	input_file=/path/to/input/PPI/file/example/ppi-$i.ppi
	./DSDmain.py --outformat matrix -o $output_file $input_file
done

# === Instructions ===
# Run this script from the directory containing DSD code.
# Specify file numbers in the first line of the for loop.
#     * Replace START and STOP with the first and last file number to execute.
#     * START and STOP are inclusive. 
# Specify the input file.
#     * Include path to the input file and the filename with extension.
# Specify the output file.
#     * Include path to the output file and the filename, WITHOUT extension.
#     * ".DSD1" or ".DSD2", etc. is automatically appended by the DSD code. 

# === Potential problems ===
# Number padding
#    Depending on the number of networks generated, filenames may appear with
#    one too many/few 0's in this script. To fix this, change the number 
#    formatting portion of the first line of the for loop. For example,
#    "%04g" will produce the string 0074 for the raw int 74. 
