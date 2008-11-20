#!/bin/sh

for i in ../../data/mutagenesis/common/s*.csv; do
	./hrmain --mutagenesis.in-file $i --out `basename $i .csv`.xml
done
