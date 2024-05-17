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
#    Copyright @ The Regents of the University of California. All rights reserved.
#

#   buildSpacedDB.sh: To build databases of discriminative spaced k-mers 
#			from the the database of discriminative 31-mers.
#


RPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

if [ ! -s "$LDIR/.dbAddress" ]; then
echo "Please run the script set_targets.sh to define the targets."
exit
fi

DIR=""
for db in `cat $LDIR/.dbAddress`
do
DIR="$db"
done

if [ ! -s $DIR/db_central_k31_t*_s1610612741_m0.tsk.sz ]; then
echo "Failed to find the database of discriminative 31-mers."
exit
fi

if [ -s $DIR/T295/db_central_k31_t*_s1610612741_m0_w22.tsk.sz ]; then
echo "Database for the first spaced seed (code name:T295) already exists."
else

$LDIR/exe/converter $DIR/db_central_k31*m0.tsk.sz 22 31 T295
if [ ! -d $DIR/T295 ]; then 
	mkdir $DIR/T295/
fi
mv $DIR/db_central_k31_t*_s1610612741_m0_w22.tsk.* $DIR/T295/
fi


if [ -s $DIR/T58570/db_central_k31_t*_s1610612741_m0_w22.tsk.sz ]; then
echo "Database for the second spaced seed (code name:T58570) already exists."
else
$LDIR/exe/converter $DIR/db_central_k31*m0.tsk.sz 22 31 T58570
if [ ! -d $DIR/T58570 ]; then
        mkdir $DIR/T58570/
fi
mv $DIR/db_central_k31_t*_s1610612741_m0_w22.tsk.* $DIR/T58570/
fi

if [ -s $DIR/T38570/db_central_k31_t*_s1610612741_m0_w22.tsk.sz ]; then
echo "Database for the third spaced seed (code name:T38570) already exists."
else
$LDIR/exe/converter $DIR/db_central_k31*m0.tsk.sz 22 31 T38570
if [ ! -d $DIR/T38570 ]; then
        mkdir $DIR/T38570/
fi
mv $DIR/db_central_k31_t*_s1610612741_m0_w22.tsk.* $DIR/T38570/
fi


# ####################################################################################
