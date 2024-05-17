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
#   Copyright 2013-2019 Ounit <rouni001@cs.ucr.edu>
#   extractSequences: To extract sequences from the input data that mapped 
#		      to a specified taxon.
#

FSCRPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

if [ $# -lt 3 ]; then
$LDIR/exe/extractSeqs
exit
fi

$LDIR/exe/extractSeqs $@

