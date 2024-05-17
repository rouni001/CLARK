/*
 * CLARK, CLAssifier based on Reduced K-mers.
 */

/*
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Copyright 2013-2019, Rachid Ounit <clark.ucr.help at gmail.com>
 */

/*
 * @author: Rachid Ounit, Ph.D Candidate.
 * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#ifndef KMERSCONVERSION_HH
#define KMERSCONVERSION_HH

#include<string>
#include<stdint.h>

//      FUNCTION
//      Name: getReverseComplement
//      Implementation notes: Convert a k-mer number into its reverse k-mer into (range of 0 to 4^k-1).
//
void getReverseComplement(const uint64_t& _ikmer, const size_t& _kmerSize, uint64_t& _rev_kmerIndex);

void getReverseComplement(const std::string& _kmer, uint64_t& _rev_kmerIndex);

//      FUNCTION
//      Name: IndexTovector
//      Implementation notes: Convert a integer into a substring (or a sequence of bases) (from 0 to 4^k-1, where d is the block size).
//
void IndexTovector(const uint64_t& _index, const size_t& _kmerSize, std::string& _kmer);

//      FUNCTION
//      Name: vectorToIndex
//      Implementation notes: Convert a substring (or sequence of bases) into an integer (range of 0 to 4^k-1).
//
void vectorToIndex(const std::string& _kmer, uint64_t& _index);

//	FUNCTION
//	Name: getSpacedSeed
//	Implementation notes: Build the spaced seed (forward) given a weight from a contiguous seed of length 31.
//
void getSpacedSeed(const std::string& _setting, const uint64_t& _kmerR, uint64_t& _skmerR);

void getSpacedSeedT38570(const uint64_t& _kmerR, uint64_t& _skmerR);
void getSpacedSeedT58570(const uint64_t& _kmerR, uint64_t& _skmerR);
void getSpacedSeedOPTSS95s2(const uint64_t& _kmerR, uint64_t& _skmerR);

#endif //KMERSCONVERSION_HH
