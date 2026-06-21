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
#  getTargetsKmers_distribution.sh: To print the target-specific k-mer distribution
#			  for the current working database. The minimum k-mer
#			  frequency defaults to 0 when omitted.
# 

LDIR=${CLARK_HOME:-$(CDPATH= cd "$(dirname "$0")/.." && pwd -P)}

if [ $# -lt 1 ]; then
echo "Usage: $0 <k-mer length: integer between 2 and 32> [min k-mers frequency: default is 0]"
exit
fi
if [ $# -eq 1 ]; then
set -- "$1" 0
fi
"$LDIR/exe/getTargetSpecificKmersStat" "$LDIR/.settings" "$@"
