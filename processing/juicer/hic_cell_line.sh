#!/usr/bin/env bash

echo "Bash version - $BASH_VERSION"

if [ $1 == "-h" ] || [ ! $1 ];
then
    echo "Script for automation of exporting HiC data using Juicebox_tools. Given a cell line name, this scripts exports all intra/inter chromosome interaction matrices.

    On usage of juicebox_tools.jar refer to http://aidenlab.org/commandlinetools/docs.html

    Usage: [-h] cell_line_name data normalization resolution bin_size output juicebox_tools_path

        * -h : show this help message and exit
        * (1) cell_line_name: string representing a cell line to work with. Available cell lines: ${!amazon_sources[@]}
        * (2) normalization: type of normalization, to be applied to the exported data. Available options: [NONE/VC/VC_SQRT/KR].
        * (3) bin_size: a size of the matrix entry in the resolution unit specified above.
        * (4) output
        * (5) juicebox_tools_path: a path to juicebox_tools.jar file
    "
    exit 0
fi

HUMAN_CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
MOUSE_CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X")
DOG_CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "X")
RABBIT_CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "X")
MACAQUE_CHROMOSOMES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "X")

declare -A cell_lines_chromosomes=(
    ["rao14_gm12878"]=HUMAN_CHROMOSOMES
    ["rao14_imr90"]=HUMAN_CHROMOSOMES
    ["rao14_hmec"]=HUMAN_CHROMOSOMES
    ["rao14_nhek"]=HUMAN_CHROMOSOMES
    ["rao14_k562"]=HUMAN_CHROMOSOMES
    ["rao14_kbm7"]=HUMAN_CHROMOSOMES
    ["rao14_huvec"]=HUMAN_CHROMOSOMES
    ["rao14_huvec"]=HUMAN_CHROMOSOMES
    ["rao14_chx12"]=MOUSE_CHROMOSOMES
    ["rudan15_dog1"]=DOG_CHROMOSOMES
    ["rudan15_dog2"]=DOG_CHROMOSOMES
    ["rudan15_mouse1"]=MOUSE_CHROMOSOMES
    ["rudan15_mouse2"]=MOUSE_CHROMOSOMES
    ["rudan15_macaque1"]=MACAQUE_CHROMOSOMES
    ["rudan15_macaque2"]=MACAQUE_CHROMOSOMES
    ["rudan15_rabbit1"]=RABBIT_CHROMOSOMES
    ["rudan15_rabbit2"]=RABBIT_CHROMOSOMES
)

for chr1 in ${cell_lines_chromosomes[$1]}; do
    for chr2 in ${cell_lines_chromosomes[$1]}; do
        ./export_hic $1 ${chr1} ${chr2} observed $2 BP $3
    done
done
