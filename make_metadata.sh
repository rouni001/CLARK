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
#   make_metadata.sh: To create meta-data for the selected database (Bacteria,
# 		      Viruses, Human or Custom). 
#

if [ $# -lt 2 ]; then

echo "Usage: $0 <Database name: bacteria, viruses, plasmid, plastid, protozoa, fungi or human> <Directory name for the database>"

exit
fi

DB=$1
DBDR=$2
TAXDR="taxonomy"
FSCRPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

if [ ! -d $DBDR ]; then
	echo "Selected directory not found. The program will create it."
	mkdir -m 775 $DBDR
	mkdir -m 775 $DBDR/Custom

	if [ ! -d $DBDR ]; then
		echo "Failed to find the directory (please check the name of directory $DBDR: Does it exist?). The program will abort."
		exit
	fi
else
	if [ ! -d $DBDR/Custom ]; then
		mkdir -m 775 $DBDR/Custom
	fi
fi

if [ ! -d $DBDR/$TAXDR ]; then
	echo "Taxonomy data missing. The program will download data to $DBDR/$TAXDR."
	mkdir -m 775 $DBDR/$TAXDR
	$LDIR/download_taxondata.sh $DBDR/$TAXDR
fi
if [ ! -f $DBDR/.taxondata ]; then
        echo "Failed to find taxonomy files. The program will try to download them..."
        $LDIR/download_taxondata.sh $DBDR/$TAXDR
        if [ ! -f $DBDR/.taxondata ]; then
                echo "Failed to find taxonomy files."
                echo "The program must abort."
                exit
        fi
fi 

if [ "$DB" != "protozoa" ] && [ "$DB" != "plastid" ] && [ "$DB" != "plasmid" ] && [ "$DB" != "fungi" ] && [ "$DB" != "custom" ] && [ "$DB" != "bacteria" ] && [ "$DB" != "viruses" ] && [ "$DB" != "human" ]; then
	echo "Failed to recognize the database: '$DB'. "
	echo "-> Supported databases are 'bacteria', 'viruses', 'plasmid', 'plastid', 'protozoa', 'fungi' and 'human' from NCBI (Please check any typo)."
	echo " Note that the choice 'bacteria' include the archaea genomes."
	echo "-> Eventually, the program can build a database from your own sequences (*.fna files) in your disk. "
	echo "If you want CLARK to use a customized database then please do the following directions: "
	echo "1) Move your sequences in fasta format with the Accession number to $DBDR/Custom/"
	echo "2) Run again this command with the option 'custom' "
	exit
fi

if [ ! -d $DBDR ]; then
	echo "The directory $DBDR does not exit. The program will create it."
	mkdir -m 775 $DBDR
	if [ ! -d $DBDR ]; then
	echo  "Failed to create the directory $DBDR. The program must abort." 
	fi
fi

if [ ! -s $DBDR/.$DB ]; then
	if [ "$DB" != "custom" ]; then
		echo "Sequences for $DB not found. The program will download them."
		$LDIR/download_RefSeqDB.sh $DBDR $DB
	else
		ls $DBDR/Custom/ > $DBDR/.$DB
		if [ ! -s $DBDR/.$DB ]; then
			echo "The database directory 'Custom' is empty."
		        echo "If you want CLARK to use a customized database then please do the following directions: "
		        echo "1) Move your sequences in fasta format with the Accession number to $DBDR/Custom/"
		        echo "2) Run again this command with the option 'custom' "
			exit
		fi
		find $DBDR/Custom/ -name '*.f*' > $DBDR/.$DB
	fi
fi

if [ ! -f $LDIR/exe/getfilesToTaxNodes ] || [ ! -f $LDIR/exe/getAccssnTaxID ]; then
	echo "Something wrong occurred (source code may be missing or unusable. Did the installation finish properly?). The program must abort."
	exit
fi

if [ ! -s $DBDR/.$DB ]; then
	echo "Failed to find the downloaded $DB sequences."
	echo "The program must abort."
	exit
fi

if [ $DB = "human" ]; then
	cd $DBDR
	if [ ! -s .$DB ]; then
		ls ./Human/*.* > .$DB
	fi
	if [ ! -s .$DB.fileToTaxIDs ]; then
		for file in `cat .$DB`
		do
		echo "$file X 9606 9605 9604 9443 40674 7711" >> .$DB.fileToTaxIDs
		done
	fi
	exit
fi

if [ ! -s $DBDR/.$DB.fileToAccssnTaxID ] ; then
	echo "Re-building $DB.fileToAccssnTaxID"
	$LDIR/exe/getAccssnTaxID $DBDR/.$DB $DBDR/$TAXDR/nucl_accss $DBDR/$TAXDR/merged.dmp > $DBDR/.$DB.fileToAccssnTaxID
fi
if [ ! -s $DBDR/.$DB.fileToTaxIDs ]; then
	echo "$DB: Retrieving taxonomy nodes for each sequence based on taxon ID..."
	$LDIR/exe/getfilesToTaxNodes $DBDR/$TAXDR/nodes.dmp $DBDR/.$DB.fileToAccssnTaxID > $DBDR/.$DB.fileToTaxIDs
fi
exit

