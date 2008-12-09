#!/bin/sh

#for i in ../../data/mutagenesis/common/s*.csv; do
	#./hrmain -c mutagen.dat -f PlainPrint --mutagenesis.in-file $i --output `basename $i .csv`.dat
#done
for i in ../../data/mutagenesis/common/s*.csv; do
	./hrmain -c mutagen.dat -f RealXMLPrint --mutagenesis.in-file $i --output `basename $i .csv`.xml
done

cp mutagen.dat ../../data/mutagenesis/results
