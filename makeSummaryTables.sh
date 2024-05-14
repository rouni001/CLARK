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
#   Copyright @ The Regents of the University of California. All rights reserved.
#
#
#
#  @author: Rachid Ounit, Ph.D.
#  @project: CLARK, Metagenomic and Genomic Sequences Classification project.
#  @note: C++/Shell IMPLEMENTATION supported on latest Linux and Mac OS.
#  makeSummaryTable.sh: To get summary tables indicating, the number of reads,
#			the number of assigned reads, and the top organisms
#			given for each CLARK report file (produced by 
#			estimate_abundance.sh) passed in parameters.

FSCRPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

if [ $# -lt 1 ]; then
echo -n "Usage: "
$LDIR/exe/makeSummaryTables 
exit
fi
$LDIR/exe/makeSummaryTables $@
