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
#   Copyright 2013-2019, Rachid Ounit <rouni001@cs.ucr.edu>
#
#   download_RefseqDB.sh: To download complete reference genomes of NCBI/RefSeq 
#   (i.e., Bacteria/Archaea, Viruses, Plasmid, Protozoa, Fungi and Human).

if [ $# -ne 2 ]; then
echo "Usage: $0 <Directory for the sequences> <Database: bacteria, viruses, plasmid, plastid, protozoa, fungi or human> "
exit
fi

FSCRPT=$(readlink -f "$0")
DIR=$(dirname "$FSCRPT")

if [ "$2" = "bacteria" ]; then

	if [ ! -s $1/.bacteria ]; then
		rm -Rf $1/Bacteria $1/.bacteria.*
		mkdir -m 775 $1/Bacteria
		cd $1/Bacteria/
		wget ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
                awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print $20}' assembly_summary.txt > .$2.tmp
		$DIR/exe/dscriptMaker .$2.tmp > ./download.sh
                chmod 711 ./download.sh
		echo "Downloading now Bacteria/Archaea complete genomes (quiet mode)... [This operation will take several hours or more to complete.]"
		./download.sh 2> .$2.tmp
		rm -f ./assembly_summary.txt
	
		wget ftp://ftp.ncbi.nih.gov/genomes/refseq/archaea/assembly_summary.txt > .$2.tmp
                awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print $20}' assembly_summary.txt > .$2.tmp 
		$DIR/exe/dscriptMaker .$2.tmp > ./download.sh
                ./download.sh 2> .$2.tmp
                rm -f .$2.tmp ./download.sh ./assembly_summary.txt

		echo "Downloading done! Uncompressing files... "
		gunzip ./*fna.gz
		find `pwd` -name '*.fna' > ../.bacteria
		cd ..
		if  [ ! -s .bacteria ]; then
			echo "Error: Failed to download bacteria sequences. "
			exit
		fi
		echo "Bacteria/Archaea sequences downloaded!"
	else
		echo "Bacteria/Archaea sequences already in $1."
	fi
	exit
fi

if [ "$2" = "viruses" ]; then
	if [ ! -s $1/.viruses ]; then
		rm -Rf $1/Viruses  $1/.viruses.*
		mkdir -m 775 $1/Viruses
		cd $1/Viruses/
		echo "Downloading now Viruses complete genomes:"
                wget ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt
                awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print $20}' assembly_summary.txt > .$2.tmp
		$DIR/exe/dscriptMaker .$2.tmp > ./download.sh
                chmod 711 ./download.sh
                ./download.sh 2> .$2.tmp
                echo "Downloading done. Uncompressing files... "
                gunzip ./*fna.gz
		rm -f .$2.tmp ./download.sh ./assembly_summary.txt

	     	find `pwd` -name '*.fna'  > ../.viruses
	     	cd ..
	      	if  [ ! -s .viruses ]; then
		      echo "Error: Failed to download viruses sequences. "
		      exit
	      	fi
	      	echo "Viruses sequences downloaded!"
	else
	      	echo "Viruses sequences already in $1."
	fi
	exit
fi

if [ "$2" = "plasmid" ]; then
        if [ ! -s $1/.plasmid ]; then
                rm -Rf $1/Plasmid  $1/.plasmid.*
                mkdir -m 775 $1/Plasmid
                cd $1/Plasmid/
                echo "Downloading now Plasmid genomes:"
		wget ftp://ftp.ncbi.nih.gov/genomes/refseq/plasmid/plasmid.1.1.genomic.fna.gz
		wget ftp://ftp.ncbi.nih.gov/genomes/refseq/plasmid/plasmid.2.1.genomic.fna.gz
		wget ftp://ftp.ncbi.nih.gov/genomes/refseq/plasmid/plasmid.3.1.genomic.fna.gz
		wget ftp://ftp.ncbi.nih.gov/genomes/refseq/plasmid/plasmid.4.1.genomic.fna.gz
                
		echo "Downloading done. Uncompressing files... "
                gunzip plasmid*.genomic.fna.gz
		echo " * Processing sequences..."
		for fileVir in `ls ./*fna`
                do
                	$DIR/exe/exeSeq $fileVir ./
                done
		rm -f ./plasmid*.genomic.fna.gz
		find `pwd` -name '*.fa'  > ../.plasmid
              	cd ..
              	if  [ ! -s .plasmid ]; then
                      echo "Error: Failed to download plasmid sequences. "
                      exit
              	fi
              	echo "Plasmid sequences downloaded!"
        else
              	echo "Plasmid sequences already in $1."
        fi
        exit
fi

if [ "$2" = "plastid" ]; then
        if [ ! -s $1/.plastid ]; then
                rm -Rf $1/Plastid  $1/.plastid.*
                mkdir -m 775 $1/Plastid
                cd $1/Plastid/
                echo "Downloading now Plastid genomes:"
                wget ftp://ftp.ncbi.nih.gov/genomes/refseq/plastid/plastid.1.1.genomic.fna.gz
                wget ftp://ftp.ncbi.nih.gov/genomes/refseq/plastid/plastid.2.1.genomic.fna.gz
                echo "Downloading done. Uncompressing files... "
		gunzip plastid*.genomic.fna.gz
		 echo " * Processing sequences..."
                for fileVir in `ls ./*fna`
                do
                        $DIR/exe/exeSeq $fileVir ./
                done
                rm -f ./plastid*.genomic.fna.gz
                find `pwd` -name '*.fa'  > ../.plastid
                cd ..
                if  [ ! -s .plastid ]; then
                      echo "Error: Failed to download plastid sequences. "
                      exit
                fi
                echo "Plastid sequences downloaded!"
        else
                echo "Plastid sequences already in $1."
        fi
        exit
fi

if [ "$2" = "fungi" ]; then
	if [ ! -s $1/.fungi ]; then
		rm -Rf $1/Fungi  $1/.fungi.*
		mkdir -m 775 $1/Fungi
		cd $1/Fungi/
		echo "Downloading now RefSeq fungi complete genomes:"
		wget ftp://ftp.ncbi.nih.gov/genomes/refseq/$2/assembly_summary.txt
		awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print $20}' assembly_summary.txt > .$2.tmp
		$DIR/exe/dscriptMaker .$2.tmp > ./download.sh
		chmod 711 ./download.sh
		./download.sh  2> .$2.tmp
		rm -f .$2.tmp ./download.sh ./assembly_summary.txt

        echo "Downloading now a list of 47 fungi reference genomes..."
#Blastomyces dermatitidis
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/525/GCA_000003525.2_BD_ER3_V1/GCA_000003525.2_BD_ER3_V1_genomic.fna.gz
#Candida parapsilosis strain CDC317
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/765/GCA_000182765.2_ASM18276v2/GCA_000182765.2_ASM18276v2_genomic.fna.gz
#Saccharomyces sp. 'boulardii'
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/413/975/GCA_001413975.1_ASM141397v1/GCA_001413975.1_ASM141397v1_genomic.fna.gz
#Blastomyces percursus strain EI222
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/883/805/GCA_001883805.1_Blas_perc_EI222_V1/GCA_001883805.1_Blas_perc_EI222_V1_genomic.fna.gz
#Aspergillus tubingensis CBS 134.48
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/890/745/GCA_001890745.1_Asptu1/GCA_001890745.1_Asptu1_genomic.fna.gz
#Fusarium solani strain JS-169
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/215/905/GCA_002215905.1_ASM221590v1/GCA_002215905.1_ASM221590v1_genomic.fna.gz
#Lomentospora prolificans strain JHH-5317
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/276/285/GCA_002276285.1_Lprolificans_pilon/GCA_002276285.1_Lprolificans_pilon_genomic.fna.gz
#Candida glabrata CBS 138
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/545/GCF_000002545.3_ASM254v2/GCF_000002545.3_ASM254v2_genomic.fna.gz
#Aspergillus fumigatus Af293
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz
#Aspergillus clavatus NRRL 1
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/715/GCF_000002715.2_ASM271v1/GCF_000002715.2_ASM271v1_genomic.fna.gz
#Ajellomyces dermatitidis SLH14081
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/855/GCF_000003855.2_ASM385v2/GCF_000003855.2_ASM385v2_genomic.fna.gz
#Aspergillus flavus NRRL3357
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/275/GCF_000006275.2_JCVI-afl1-v2.0/GCF_000006275.2_JCVI-afl1-v2.0_genomic.fna.gz
#Candida tropicalis MYA-3404
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/335/GCF_000006335.2_ASM633v2/GCF_000006335.2_ASM633v2_genomic.fna.gz
#Candida dubliniensis CD36
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/945/GCF_000026945.1_ASM2694v1/GCF_000026945.1_ASM2694v1_genomic.fna.gz
#Cryptococcus neoformans var. grubii H99
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz
#Coccidioides immitis RS
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/335/GCF_000149335.2_ASM14933v2/GCF_000149335.2_ASM14933v2_genomic.fna.gz
#Meyerozyma guilliermondii ATCC 6260
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/425/GCF_000149425.1_ASM14942v1/GCF_000149425.1_ASM14942v1_genomic.fna.gz
#Fusarium verticillioides 7600
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/555/GCF_000149555.1_ASM14955v1/GCF_000149555.1_ASM14955v1_genomic.fna.gz
#Histoplasma capsulatum NAm1
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/585/GCF_000149585.1_ASM14958v1/GCF_000149585.1_ASM14958v1_genomic.fna.gz
#Fusarium oxysporum f. sp. lycopersici 4287
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/955/GCF_000149955.1_ASM14995v2/GCF_000149955.1_ASM14995v2_genomic.fna.gz
#Paracoccidioides lutzii Pb01
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/705/GCF_000150705.2_Paracocci_br_Pb01_V2/GCF_000150705.2_Paracocci_br_Pb01_V2_genomic.fna.gz
#Paracoccidioides brasiliensis Pb18
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/735/GCF_000150735.1_Paracocci_br_Pb18_V2/GCF_000150735.1_Paracocci_br_Pb18_V2_genomic.fna.gz
#Coccidioides posadasii C735 delta SOWgp
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/335/GCF_000151335.1_JCVI-cpa1-1.0/GCF_000151335.1_JCVI-cpa1-1.0_genomic.fna.gz
#Candida albicans SC5314
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz
#Pneumocystis murina B123
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/349/005/GCF_000349005.2_Pneumo_murina_B123_V4/GCF_000349005.2_Pneumo_murina_B123_V4_genomic.fna.gz
#[Candida] auris strain 6684
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/189/475/GCF_001189475.1_ASM118947v1/GCF_001189475.1_ASM118947v1_genomic.fna.gz
#Pneumocystis jirovecii RU7
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/477/535/GCF_001477535.1_Pneu_jiro_RU7_V2/GCF_001477535.1_Pneu_jiro_RU7_V2_genomic.fna.gz
#Pneumocystis carinii B80
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/477/545/GCF_001477545.1_Pneu_cari_B80_V3/GCF_001477545.1_Pneu_cari_B80_V3_genomic.fna.gz
#Pichia kudriavzevii strain 129
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/983/325/GCF_001983325.1_ASM198332v1/GCF_001983325.1_ASM198332v1_genomic.fna.gz
#Aspergillus fischeri
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/645/GCF_000149645.1_ASM14964v1/GCF_000149645.1_ASM14964v1_genomic.fna.gz
#Fusarium proliferatum
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/067/095/GCA_900067095.1_F._proliferatum_ET1_version_1/GCA_900067095.1_F._proliferatum_ET1_version_1_genomic.fna.gz
#Malassezia spp.
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/264/985/GCA_001264985.1_ASM126498v1/GCA_001264985.1_ASM126498v1_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/600/835/GCA_001600835.1_JCM_12085_assembly_v001/GCA_001600835.1_JCM_12085_assembly_v001_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/600/795/GCA_001600795.1_JCM_11963_assembly_v001/GCA_001600795.1_JCM_11963_assembly_v001_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/600/775/GCA_001600775.1_JCM_11348_assembly_v001/GCA_001600775.1_JCM_11348_assembly_v001_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/551/515/GCA_002551515.1_Malafurf/GCA_002551515.1_Malafurf_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/264/725/GCA_001264725.1_ASM126472v1/GCA_001264725.1_ASM126472v1_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/181/695/GCF_000181695.1_ASM18169v1/GCF_000181695.1_ASM18169v1_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/349/305/GCF_000349305.1_ASM34930v2/GCF_000349305.1_ASM34930v2_genomic.fna.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/264/965/GCA_001264965.1_ASM126496v1/GCA_001264965.1_ASM126496v1_genomic.fna.gz
#Cladosporium sphaerospermum
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/261/425/GCA_000261425.2_UM843Vel2.0/GCA_000261425.2_UM843Vel2.0_genomic.fna.gz
#Candida orthopsilosis
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/315/875/GCF_000315875.1_ASM31587v1/GCF_000315875.1_ASM31587v1_genomic.fna.gz
#Alternaria alternata
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/642/055/GCF_001642055.1_Altal1/GCF_001642055.1_Altal1_genomic.fna.gz
#Cryptococcus neoformans
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz
#Trichosporon asahii
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/293/215/GCF_000293215.1_Trichosporon_asahii_1/GCF_000293215.1_Trichosporon_asahii_1_genomic.fna.gz
#Enterocytozoon bieneusi
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/485/GCF_000209485.1_ASM20948v1/GCF_000209485.1_ASM20948v1_genomic.fna.gz
#Naganishia albida
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/599/735/GCA_001599735.1_JCM_2334_assembly_v001/GCA_001599735.1_JCM_2334_assembly_v001_genomic.fna.gz
		echo "Downloading done. Uncompressing files... "
		gunzip ./*fna.gz

		find `pwd` -name '*.fna' > ../.fungi
	  	cd ../
	  	if  [ ! -s .fungi ]; then
	  		echo "Error: Failed to download fungi sequences. "
	  		exit
	  	fi
	  	echo "Fungi sequences downloaded!"
	  else
	  	echo "Fungi sequences already in $1."
	  fi
	  exit
fi

if [ "$2" = "protozoa" ]; then
        if [ ! -s $1/.protozoa ]; then
                rm -Rf $1/Protozoa  $1/.protozoa.*
                mkdir -m 775 $1/Protozoa
                cd $1/Protozoa/
                echo "Downloading now RefSeq protozoa complete genomes:"
                wget ftp://ftp.ncbi.nih.gov/genomes/refseq/$2/assembly_summary.txt
                awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print $20}' assembly_summary.txt > .$2.tmp
                $DIR/exe/dscriptMaker .$2.tmp > ./download.sh
                chmod 711 ./download.sh
                ./download.sh 2> .$2.tmp
                rm -f .$2.tmp ./download.sh ./assembly_summary.txt
                echo "Downloading done. Uncompressing files... "
                gunzip ./*fna.gz

                find `pwd` -name '*.fna' > ../.protozoa
                cd ../
                if  [ ! -s .protozoa ]; then
                        echo "Error: Failed to download protozoa sequences. "
                        exit
                fi
                echo "Protozoa sequences downloaded!"
          else
                echo "Protozoa sequences already in $1."
          fi
          exit
fi

if [ "$2" = "human" ]; then
	  if [ ! -s $1/.human ]; then
	  	rm -Rf $1/Human  $1/.human.*
	  	mkdir -m 775 $1/Human
	  	cd $1/Human/
	  	echo "Downloading now latest Human genome:"
	  	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.37_GRCh38.p11/GCF_000001405.37_GRCh38.p11_genomic.fna.gz
		echo "Downloading done. Uncompressing files... "
	  	gunzip ./*fna.gz

	  	find `pwd` -name '*.fna' > ../.human
	  	cd ../
	  	if  [ ! -s .human ]; then
	  		echo "Error: Failed to download human sequences. "
	  		exit
	  	fi
	  	echo "Human genome downloaded!"
	else
		echo "Human genome already in $1."
	fi
	exit
fi

echo "Failed to recognize parameter: $2. Please choose between: bacteria, viruses, plasmid, plastid, protozoa, fungi or human."

