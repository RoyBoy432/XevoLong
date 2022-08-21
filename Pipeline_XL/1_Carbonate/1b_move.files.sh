#!/bin/bash

cd ~/XevoLong.lnk

AR=( $(seq 1 88 ) )

for i in "${AR[@]}"

do
	cp *XL${i}_S* ./Sample_${i}
	
	done