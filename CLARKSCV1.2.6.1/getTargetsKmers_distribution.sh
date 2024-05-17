#! /bin/sh

#
# CLARK, CLAssifier based on Reduced K-mers.
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
#
#  @author: Rachid Ounit, Ph.D Candidate.
#  @project: CLARK, Metagenomic and Genomic Sequences Classification project.
#  @note: C++/Shell IMPLEMENTATION supported on latest Linux and Mac OS.
#  getTargetsKmers_distribution.sh: To estimate abundance of target identified (reported
#			  by taxa name and lineage) with count and proportion 
#			  (against all reads or classified reads). 
#			  Filtering options are offered.
# 

FSCRPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

if [ $# -lt 2 ]; then
echo "Usage: $0 <k-mer length: integer between 2 and 32> <min k-mers frequency: default is 0>"
exit
fi
$LDIR/exe/getTargetSpecificKmersStat $LDIR/.settings $@

