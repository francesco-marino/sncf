#!/bin/bash

int_name="";
flag=0
c0=0; c1=0
n_args=$#
if [[ $n_args -ge 1 ]]; then
   int_name=$1
   flag=1
   if [[ $n_args -ge 3 ]]; then
       c0=$2; c1=$3
       echo -e "Gradient terms introduced: C0=$c0\tC1=$c1"
   fi
else
   echo "Reading file sncf.in"
#  read 2nd line of sncf.in 
   int_name=$(sed -n '2{p;q;}' sncf.in)
fi

int_name=${int_name,,}
if [[ $flag -ge 1 ]]; then
  file="Interactions/${int_name}.in"
  if [ -f "$file" ]; then
    cp "${file}" "interaction.in"
  else
    echo "Error: interaction not available"
    exit 1
  fi
fi


if [[ $c0 -ne 0 || $c1 -ne 0 ]]; then
   int_name="${int_name}_${c0#-}_${c1#-}"
fi
  

nuclei=(o si s ca ca ca ca zr sn pb)
A=(16 34 36 36 40 48 52 90 132 208)
Z=(8  14 16 20 20 20 20 40 50  82)
size=${#A[*]}

echo -e  "Using interaction ${int_name}\n"
# write nucl, en, radius
tab_file="Results/tab_${int_name}.dat"


for ((i=0; i<size;i++)); do
#  name="${nuclei[$i]}${A[$i]}"
   name="${A[$i]}${nuclei[$i]}"
   echo -e "\nStudying nucleus ${name}"
   ./sncf.x ${A[$i]} ${Z[$i]} ${c0} ${c1} | tee "Results/summary_${name}_${int_name}.dat"
   
   cp "dens.out" "Results/rho_${name}_${int_name}_hf.dat"
   cp "energy_density_full.out" "Results/en_dens_${name}_${int_name}_hf.dat" 

   if [[ $flag -eq 0 ]]; then
      cp "contributions.out" "Results/contr_${name}_${int_name}.dat"
   fi

   t=$(<"temp.out")
   a=($(echo "$t" | tr ' ' '\n'))
   if [[ $i -eq 0 ]] 
   then
     echo -e "#Nucl\tEnergy\t\tR(ch)" > ${tab_file}
   fi
   echo -e "${name}\t${a[2]}\t${a[3]}" >> ${tab_file}
done

echo -e "\nDone! Look inside the Results folder"

