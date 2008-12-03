#!/bin/sh

for i in ../../data/mutagenesis/common/s*.csv; do
	./hrmain -c mutagen.dat --mutagenesis.in-file $i --output `basename $i .csv`.xml
done
