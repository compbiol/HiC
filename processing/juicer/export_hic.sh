#!/usr/bin/env bash

echo "Bash version - $BASH_VERSION"

declare -A amazon_sources=(
    ["rao14_gm12878"]="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic"
    ["rao14_imr90"]="https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic"
    ["rao14_hmec"]="https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic"
    ["rao14_nhek"]="https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic"
    ["rao14_k562"]="https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic"
    ["rao14_kbm7"]="https://hicfiles.s3.amazonaws.com/hiseq/kbm7/in-situ/combined.hic"
    ["rao14_huvec"]="https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic"
    ["rao14_hela"]="https://hicfiles.s3.amazonaws.com/hiseq/hela/in-situ/combined.hic"
    ["rao14_chx12-lx"]="https://hicfiles.s3.amazonaws.com/hiseq/ch12-lx-b-lymphoblasts/in-situ/combined.hic"
    ["rudan15_dog1"]="http://hicfiles.s3.amazonaws.com/external/rudan/canis-lupus-rep1.hic"
    ["rudan15_dog2"]="http://hicfiles.s3.amazonaws.com/external/rudan/canis-lupus-rep1.hic"
    ["rudan15_mouse1"]="http://hicfiles.s3.amazonaws.com/external/rudan/mouse-rep1.hic"
    ["rudan15_mouse2"]="http://hicfiles.s3.amazonaws.com/external/rudan/mouse-rep2.hic"
    ["rudan15_macaque1"]="http://hicfiles.s3.amazonaws.com/external/rudan/macaque-rep1.hic"
    ["rudan15_macaque2"]="http://hicfiles.s3.amazonaws.com/external/rudan/macaque-rep2.hic"
    ["rudan15_rabbit1"]="http://hicfiles.s3.amazonaws.com/external/rudan/rabbit-rep1.hic"
    ["rudan15_rabbit2"]="http://hicfiles.s3.amazonaws.com/external/rudan/rabbit-rep2.hic"
)

if [ $1 == "-h" ] || [ ! $1 ];
then
    echo "Script for automation of exporting HiC data using Juicebox_tools

    On usage of juicebox_tools.jar refer to http://aidenlab.org/commandlinetools/docs.html

    Usage: [-h] cell_line_name chr1 chr2 data normalization resolution bin_size output juicebox_tools_path

        * -h : show this help message and exit
        * (1) cell_line_name: string representing a cell line to work with. Available cell lines: ${!amazon_sources[@]}
        * (2) chr1: name of the first chromosome in the exported matrix dataset
        * (3) chr2: name of the second chromosome in the exported matrix dataset
        * (4) data: a type of data, one wishes to export. Available options: [observed|oe|pearson|norm|expected], where last two produce vectors.
        * (5) normalization: type of normalization, to be applied to the exported data. Available options: [NONE/VC/VC_SQRT/KR].
        * (6) resolution: a unit of resolution. Available options: [BP|FRAG]
        * (7) bin_size: a size of the matrix entry in the resolution unit specified above.
        * (8) juicebox_tools_path: a path to juicebox_tools.jar file
    "
    exit 0
fi


if [ ! ${amazon_sources[$1]+_} ];
then
     echo "Cell line \"$1\" is not present in available sources"
     echo -e "Available cell lines: \n\t ${!amazon_sources[@]}"
     exit 1
fi

JUICEBOX_TOOLS_PATH=$8
if [[ ! $8 ]];
then
    JUICEBOX_TOOLS_PATH="juicebox_tools.jar"
fi

OUTPUT_FILE_NAME="$1_$2_$3_$7_$5.txt"
echo "Result file: $OUTPUT_FILE_NAME"
# java -jar juicebox_tools.jar dump observed file.hic chr1 chr2 BP 10000 output.txt
java -jar ${JUICEBOX_TOOLS_PATH} dump $4 $5 ${amazon_sources[$1]} $2 $3 $6 $7 ${OUTPUT_FILE_NAME}

