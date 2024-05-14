/*
 * CONVERTER (CONTIGUOUS SEED TO SPACED SEED), program part of 
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

#ifndef  CONTIGUOUSTOSPACED_HH
#define  CONTIGUOUSTOSPACED_HH

#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <cstring>

#include <stdint.h>
#include "./hashTable_hh.hh"
#include "./dataType.hh"

template <typename HKMERr, typename HKMERrs>
class contiguousTospaced
{
	private:
		hTable<HKMERrs,lElement> 	m_sTable;
		size_t				m_weight;
		size_t				m_len;
		string				m_setting;
	public:
		contiguousTospaced():m_sTable(), m_weight(22), m_len(31), m_setting("T295"){}

		contiguousTospaced(const size_t& _weight, const std::string& _setting, const size_t& _len = 31):m_sTable(), m_weight(_weight), m_setting(_setting), m_len(_len)	
	{
	}

		~contiguousTospaced(){}

		void populate(const char* _filename);
		void write(const char* _filename);
	private:
		bool update(const uint64_t& _km, const ILBL& _lbl);

};

#endif //  CONTIGUOUSTOSPACED_HH

#include "./kmersConversion.hh"

	template <typename HKMERr, typename HKMERrs>
void contiguousTospaced<HKMERr,HKMERrs>::populate(const char* _filename)
{
	char * file_lbl = (char*) calloc(strlen(_filename)+4,sizeof(char));
	char * file_key = (char*) calloc(strlen(_filename)+4,sizeof(char));
	char * file_sze = (char*) calloc(strlen(_filename)+4,sizeof(char));

	sprintf(file_lbl, "%s.lb", _filename);
	sprintf(file_key, "%s.ky", _filename);
	sprintf(file_sze, "%s.sz", _filename);
	FILE * fd_l = fopen(file_lbl,"r");
	FILE * fd_k = fopen(file_key,"r");
	FILE * fd_s = fopen(file_sze,"r");

#define LEN 100000

	if (fd_l == NULL)
	{       std::cerr << "Failed to open " << file_lbl << std::endl; return;    }
	if (fd_k == NULL)
	{       std::cerr << "Failed to open " << file_key << std::endl; return;   }
	if (fd_s == NULL)
	{       std::cerr << "Failed to open " << file_sze << std::endl; return;   }


	ITYPE t = 0, i  = 0;

	uint8_t         c[LEN];
	ILBL            lbl[LEN];
	HKMERr          key[LEN];
	size_t len      = fread(c,   1, LEN, fd_s);
	size_t len_l    = fread(lbl, sizeof(ILBL), LEN, fd_l);
	size_t len_k    = fread(key, sizeof(HKMERr), LEN, fd_k);

	size_t nbElement = 0, u = 0, v = 0, v_c = 0;

	std::cerr << "Conversion on-going... " << std::endl;

	while (true)
	{
		while (i < len)
		{
			if (c[i] > 0)
			{
				nbElement += c[i];
				u = 0;
				v_c = v;
				while (u < c[i])
				{
					if (v_c < len_l)
					{
						uint64_t km = ((uint64_t) t + ((uint64_t) HTSIZE * ((uint64_t) key[v_c])));
						ILBL albl = lbl[v_c++];
						u++;
						update(km, albl);
						continue;
					}
					len_l = fread(lbl, sizeof(ILBL), LEN, fd_l);
					len_k = fread(key, sizeof(HKMERr), LEN, fd_k);
					v_c = 0;
					if (len_l == 0)
					{       break;  }
				}
				v = v_c;
			}
			i++;
			t++;
		}
		len = fread(c, 1, LEN, fd_s);
		i = 0;
		if (len == 0)
		{       break;  }
	}
	fclose(fd_l);
	fclose(fd_k);
	fclose(fd_s);

	free(file_lbl);
	file_lbl=NULL;
	free(file_key);
	file_key=NULL;
	free(file_sze);
	file_sze=NULL;
	std::cerr << nbElement << " contiguous k-mers were successfully processed. " << std::endl;
}

	template <typename HKMERr, typename HKMERrs>
bool contiguousTospaced<HKMERr,HKMERrs>::update(const uint64_t& _km, const ILBL& _lbl)
{
	uint64_t skm = 0;
	// Try forward
	getSpacedSeed(m_setting, _km, skm);
	size_t xE, yE;
	ILBL h;
	IOCCR mult;
	ICount count;

	if (m_sTable.find(skm, xE, yE, h, mult, count))
	{
		if (_lbl != h)
		{
			m_sTable.unmarkElementAt(xE, yE);	
			return false;
		}
		return true;
	}
	// Try reverse
	uint64_t _ikmerR = _km, skmr =0;
	
	getReverseComplement(_km, m_len, _ikmerR);
	getSpacedSeed(m_setting, _ikmerR, skmr);
	if (m_sTable.find(skmr, xE, yE, h, mult, count))
	{
		if (_lbl != h)
		{
			m_sTable.unmarkElementAt(xE, yE);
			return false;
		}
		return true;
	}
	m_sTable.insertMarked(skm, _lbl, 1);
	return true;
}

	template <typename HKMERr, typename HKMERrs>
void contiguousTospaced<HKMERr,HKMERrs>::write(const char* _filename)
{
	char * newfilename = (char*) calloc(strlen(_filename)+5,sizeof(char));
	for(size_t t = 0; t < strlen(_filename) - 4; t++)
	{
		newfilename[t] = _filename[t];
	}
	sprintf(newfilename, "%s_w%lu.tsk", newfilename, m_weight);
	m_sTable.sortall(2);
	m_sTable.write(newfilename, 2, false);
}
