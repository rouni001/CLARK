#!/bin/sh

# 
#   CLARK, CLAssifier based on Reduced K-mers.
# 
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#   Copyright @ The Regents of the University of California. All rights reserved.
#
#   set_targets.sh: To create targets definition of selected databases
#            (Bacteria, Viruses, Plasmid, Protozoa, Fungi Human and Custom).

if [ $# -lt 2 ]; then

echo "Usage: $0 <Directory_path> <Database_choice+: bacteria, viruses, plasmid, plastid, protozoa, fungi, human, custom. Recommended database is: bacteria viruses> <taxonomy rank: --phylum, --class, --order, --family, --genus or --species. Default is: --species>"

exit
fi

DBDR=$1
RANK=0
FSCRPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

if [ ! -d $DBDR ]; then
	echo "Selected directory not found. The program will create it."
	mkdir -m 775 $DBDR
	if [ ! -d $DBDR ]; then
	echo "Failed to create the directory (please check the name of directory $DBDR and whether it exists). The program will abort."
	exit
	fi
fi
echo $DBDR > $LDIR/.DBDirectory
for var in $@
do
      if [ "$var" = "--species" ]; then
        RANK=0
        break
        fi
        if [ "$var" = "--genus" ]; then
        RANK=1
        break
        fi
        if [ "$var" = "--family" ]; then
        RANK=2
        break
        fi
        if [ "$var" = "--order" ]; then
        RANK=3
        break
        fi
        if [ "$var" = "--class" ]; then
        RANK=4
        break
        fi
        if [ "$var" = "--phylum" ]; then
        RANK=5
        break
        fi
	PREF=`echo $var | cut -c1-2`
        if [ "$PREF" = "--" ]; then
		echo "Failed to recognize this parameter: $var"
		exit 
	fi
done

if [ -f $DBDR/targets.txt ]; then
	rm -f $DBDR/targets.txt
fi

touch $DBDR/targets.txt
rm -f $DBDR/.tmp $LDIR/.settings $LDIR/files_excluded.txt $DBDR/files_excluded.txt
subDB=""
us="_"
for db in $@
do
	if [ "$db" != "$DBDR" ]; then
		PRE=`echo $db | cut -c1-2`
		if [ "$PRE" != "--" ]; then

			echo -n "Collecting metadata of $db... "
			$LDIR/make_metadata.sh $db $DBDR
			if [ ! -s $DBDR/.$db ]; then
				exit
			fi
			if [ ! -f $DBDR/.taxondata ]; then
				exit
			fi
			echo "done."
			if [ -s $DBDR/.$db.fileToTaxIDs ]; then 
				$LDIR/exe/getTargetsDef $DBDR/.$db.fileToTaxIDs $RANK >> $DBDR/targets.txt 
				subDB="$subDB$db$us"
				if [ -f $LDIR/files_excluded.txt ]; then
					cat $LDIR/files_excluded.txt >> $DBDR/.tmp
					rm $LDIR/files_excluded.txt
				fi
			fi
		fi
	fi
done

subDB="$subDB$RANK"
echo "-T $DBDR/targets.txt" > $LDIR/.settings
if [ ! -d $DBDR/$subDB ]; then
	echo "Creating directory to store discriminative k-mers: $DBDR/$subDB"
	mkdir -m 775 $DBDR/$subDB
fi
echo "-D $DBDR/$subDB/" >> $LDIR/.settings
if [ -s $DBDR/.tmp ]; then
	mv $DBDR/.tmp $DBDR/files_excluded.txt
fi
echo "$DBDR/$subDB" > $LDIR/.dbAddress
