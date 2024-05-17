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

   Copyright @The Regents of the University of California. All rights reserved.

*/
/*
 * @author: Rachid Ounit, Ph.D.
 * @project: CLARK, Metagenomics and Genomics Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
using namespace std;

#include "kmersConversion.hh"
#include "dataType.hh"

// Code from Jellyfish
void getReverse(uint64_t& _km_r, uint8_t k)
{
	_km_r = ((_km_r >> 2)  & 0x3333333333333333UL) | ((_km_r & 0x3333333333333333UL) << 2);
	_km_r = ((_km_r >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((_km_r & 0x0F0F0F0F0F0F0F0FUL) << 4);
	_km_r = ((_km_r >> 8)  & 0x00FF00FF00FF00FFUL) | ((_km_r & 0x00FF00FF00FF00FFUL) << 8);
	_km_r = ((_km_r >> 16) & 0x0000FFFF0000FFFFUL) | ((_km_r & 0x0000FFFF0000FFFFUL) << 16);
	_km_r = ( _km_r >> 32                        ) | ( _km_r                        << 32);
	_km_r = (((uint64_t)-1) - _km_r) >> (64 - (k << 1));
}

void getKmers(const std::string& c, uint64_t& _km_f, uint8_t k)
{
	_km_f = 0;
	uint8_t cpt = 0;
	while (cpt < k)
	{
		_km_f <<=2;
		switch (c[cpt++])
		{
			case 'A': case 'a':_km_f ^= 3; break;
			case 'C': case 'c':_km_f ^= 2; break;
			case 'G': case 'g':_km_f ^= 1; break;
			case 'T': case 't':_km_f ^= 0; break;
			default:
					   cerr << "Failed to compute k-mer value of " << c << endl;
					   std::exit(1);
					   break;
		}
	}
}

void getReverseComplement(const uint64_t& _ikmer, const size_t& _kmerSize, uint64_t& _rev_kmerIndex)
{
	_rev_kmerIndex = _ikmer;
	getReverse(_rev_kmerIndex, (uint8_t)_kmerSize);
}

void getReverseComplement(const std::string& _kmer, uint64_t& _rev_kmerIndex)
{
	_rev_kmerIndex = 0;
	getKmers(_kmer, _rev_kmerIndex,  (uint8_t)_kmer.size());
	getReverse(_rev_kmerIndex, (uint8_t)_kmer.size());
}

void vectorToIndex(const std::string& _kmer, uint64_t& _index)
{
	getKmers(_kmer, _index, _kmer.size());
}

void IndexTovector(const uint64_t& _index, const size_t& _kmerSize, std::string& _kmer)
{
	uint64_t num(_index);
	uint64_t max = 1;
	for(size_t t = 0 ; t < _kmerSize - 1; t++)
	{
		max *= 4;
	}
	vector<char> subKmer;
	while (max > 0)
	{
		if (num >= max)
		{
			if (num/max == 1)
			{
				subKmer.push_back('G');
				num = num - max;
			}
			if (num/max == 2)
			{
				subKmer.push_back('C');
				num = num - 2*max;
			}
			if (num/max == 3)
			{
				subKmer.push_back('A');
				num = num - 3*max;
			}
		}
		else
		{
			subKmer.push_back('T');
		}
		max = max / 4;
	}
	//Copy and reversing char order
	string k = "";
	for(size_t t = 0 ; t < subKmer.size() ; t++)
	{
		//k.push_back(subKmer[subKmer.size() - 1 - t]);
		k.push_back(subKmer[t]);
	} 
	_kmer = k;
}

/*
 *
 1,1111101101001110100111011101111, hitProba: 0.998113
 1,1111011101110010111001011011111, hitProba: 0.998113
 2,1111100101110110101100111011111, hitProba: 0.998099
 2,1111101110011010110111010011111, hitProba: 0.998099
 3,1111101110011011010111001011111, hitProba: 0.998093
 3,1111101001110101101100111011111, hitProba: 0.998093
 4,1111101011100101101110011011111, hitProba: 0.99809
 4,1111101100111011010011101011111, hitProba: 0.99809
 5,1111011101011011001110100111111, hitProba: 0.998084
 5,1111110010111001101101011101111, hitProba: 0.998084
 *
 */

void getSpacedSeed(const std::string& _setting, const uint64_t& _kmerR, uint64_t& _skmerR)
{
	if (_setting == "T295")
	{	return getSpacedSeedOPTSS95s2(_kmerR, _skmerR); }
	if (_setting == "T58570")
	{       return getSpacedSeedT58570(_kmerR, _skmerR); }
	if (_setting == "T38570")
	{       return getSpacedSeedT38570(_kmerR, _skmerR); }
	cerr << "Failed to find the mask for the spaced seed requested."<< endl;
	exit(1);
}

void getSpacedSeedT58570(const uint64_t& _kmerR,  uint64_t& _skmerR)
{
	// 1111101001110101101100111011111
	// 11111*1**111*1*11*11**111*11111
	_skmerR = _kmerR >> 52; // 26, Add 1
	_skmerR <<= 2 ;
	_skmerR ^= (_kmerR & 0x3FFFFFFFFFFFFUL)>> 48; // Add 1
	_skmerR <<= 6 ;
	_skmerR ^= (_kmerR & 0xFFFFFFFFFFFUL)>> 38; // Add 1
	_skmerR <<= 2 ;
	_skmerR ^= (_kmerR & 0xFFFFFFFFFUL)>> 34; // Add 1
	_skmerR <<= 4 ;
	_skmerR ^= (_kmerR & 0xFFFFFFFFUL)>> 28; // Add 1
	_skmerR <<= 4 ;
	_skmerR ^= (_kmerR & 0x3FFFFFFUL)>> 22; // Add 1
	_skmerR <<= 6 ;
	_skmerR ^= (_kmerR & 0x3FFFFUL)>> 12; // Add 1
	_skmerR <<= 10;
	_skmerR ^= (_kmerR & 0x3FFUL); // Add 1
	return;
}
//////////
void getSpacedSeedT38570(const uint64_t& _kmerR, uint64_t& _skmerR)
{
	// 1111101011100101101110011011111
	// 11111*1*111**1*11*111**11*11111
	_skmerR = _kmerR >> 52; // 26, Add 1
	_skmerR <<= 2 ;
	_skmerR ^= (_kmerR & 0x3FFFFFFFFFFFFUL)>> 48; // Add 1
	_skmerR <<= 6 ;
	_skmerR ^= (_kmerR & 0x3FFFFFFFFFFFUL)>> 40; // Add 1
	_skmerR <<= 2 ;
	_skmerR ^= (_kmerR & 0xFFFFFFFFFUL)>> 34; // Add 1
	_skmerR <<= 4 ;
	_skmerR ^= (_kmerR & 0xFFFFFFFFUL)>> 28; // Add 1
	_skmerR <<= 6 ;
	_skmerR ^= (_kmerR & 0x3FFFFFFUL)>> 20; // Add 1
	_skmerR <<= 4 ;
	_skmerR ^= (_kmerR & 0xFFFFUL)>> 12; // Add 1
	_skmerR <<= 10;
	_skmerR ^= (_kmerR & 0x3FFUL); // Add 1
	return;
}

void getSpacedSeedOPTSS95s2(const uint64_t& _kmerR,  uint64_t& _skmerR)
{
	// 1111*111*111**1*111**1*11*11111
	_skmerR = _kmerR >> 54; // Add 1111
	_skmerR <<= 6 ;
	_skmerR ^= (_kmerR & 0xFFFFFFFFFFFFFUL)>> 46; // Add 111
	_skmerR <<= 6 ;
	_skmerR ^= (_kmerR & 0xFFFFFFFFFFFUL) >> 38; // Add 111
	_skmerR <<= 2 ;
	_skmerR ^= (_kmerR & 0x3FFFFFFFFUL) >> 32; // Add 1
	_skmerR <<= 6 ;
	_skmerR ^= (_kmerR & 0x3FFFFFFFUL) >> 24; // Add 111
	_skmerR <<= 2 ;
	_skmerR ^= (_kmerR & 0xFFFFFUL) >> 18; // Add 1
	_skmerR <<= 4 ;
	_skmerR ^= (_kmerR & 0xFFFFUL) >> 12; // Add 11
	_skmerR <<= 10 ;
	_skmerR ^= (_kmerR & 0x3FFUL); // Add 11111
	return;
}

