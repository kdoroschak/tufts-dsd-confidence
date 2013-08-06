#!/bin/bash

# DsdRunner.sh
# Run a batch of files through DSD

# See below the script for instructions.

for i in $(seq -f "%03g" 1 1)
do
	output_path=~/bcb/Biogrid/biogrid_gosemsim/biogrid_ppi_gosemsim_transition
	output_name=biogrid_ppi_gosemsim_transition
	input_file=~/bcb/Biogrid/biogrid_gosemsim/biogrid_ppi_nom.gosemsim.prob
	transition_matrix=~/bcb/Biogrid/biogrid_gosemsim/biogrid_ppi_nom.gosemsim.prob.matrix
	output_file=$output_path"/"$output_name
	if [ ! -d "$output_path" ]; then
		echo "Output directory does not exist. You entered:"
		echo $output_path
		exit 1
	fi

	./DSDmain.py --outformat matrix -o $output_file -tm $transition_matrix $input_file
	if [ $? -eq 0 ]; then
		chmod 664 $output_file.DSD1
		chgrp bcb $output_file.DSD1
	fi
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
