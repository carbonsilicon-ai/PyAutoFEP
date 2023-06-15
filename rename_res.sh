#!/bin/bash
#

input=$1
sedcommands=( "s/HID /HISD/" "s/HIE /HISE/" "s/HIS /HISH/" "s/HB3/HB1/" "s/HG3/HG1/" "s/HG13 ILE/HG11 ILE/" "s/HE3 LYS/HE1 LYS/" "s/HA3 GLY/HA1 GLY/" "s/HD3 PRO/HD1 PRO/" "s/HD3 LYS/HD1 LYS/" "s/HD3 ARG/HD1 ARG/" )

[ -z "$input" ] && { printf "Rename some residues\nUsage: $0 input.pdb\n"; exit 1; }

cp "${input}" "${input}.bak"

for ((i = 0; i < ${#sedcommands[@]}; i++))
do
    command="${sedcommands[$i]}"
    sed -i "$command" "$input"
done
