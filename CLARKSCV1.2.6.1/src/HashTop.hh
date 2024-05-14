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

#ifndef HASHTOP_HH
#define HASHTOP_HH

#define MAXLEN	65536

#include "./dataType.hh"
#include<string.h>
#include<sstream>

class HashTop
{
	public:
	
	HashTop():m_ISize(0), m_Token(0), m_IBest(0)
	{
		for(size_t i =0; i < MAXLEN; i++)
		{	
			m_ITable[i] = 0;
			m_CTable[i] = 0;
		}
		init();
	}
	~HashTop()
	{
	}
	void init()
	{
		m_Token = 1;
		m_ISize = 0;
                m_IBest = 0;
		m_CBest = 0;
		m_Total = 0;
	}
	void next()
	{
		m_Token++;
		m_ISize = 0;
		m_IBest = 0;
		m_CBest = 0;
		m_Total = 0;
	}
	void insert(const ILBL& _e)
	{
		if (m_ITable[_e] == m_Token)
		{
			m_CTable[_e]++;
		}
		else
		{
			m_ITable[_e] = m_Token;
			m_CTable[_e] = 1;
			m_Indexes[m_ISize++] = _e;
		}
		if (m_CTable[_e] > m_CBest)
                {
                	m_IBest = _e;
			m_CBest = m_CTable[_e];
                }
		m_Total++;
	}
	void insert(const ILBL& _e, const size_t& _m)
        {
                if (m_ITable[_e] == m_Token)
                {
                        m_CTable[_e] += _m;
                }
                else
                {
                        m_ITable[_e] = m_Token;
                        m_CTable[_e] = _m;
                        m_Indexes[m_ISize++] = _e;
                }
                if (m_CTable[_e] > m_CBest)
                {
			m_CBest = m_CTable[_e];
                        m_IBest = _e;
                }
		m_Total += (ITYPE) _m;
        }
	void getBest(ITYPE& _h, ITYPE& _c)
	{
		if (m_ISize == 0)
		{
			_h = 0;
			_c = 0;
			return;
		}
		_h = m_IBest;
		_c = m_CBest;
	}
	void getSecondBest(ITYPE& _h, ITYPE& _c)
	{
		if (m_ISize == 0)
                {
                        _h = 0;
                        _c = 0;
			return;
                }
		_h = 0;
		ITYPE sbest = 0;
		ITYPE tmp = 0;
		for(size_t t = 0; t < m_ISize; t++)
		{
			tmp = m_Indexes[t];
			if (tmp == m_IBest)
				continue;
			if (m_CTable[tmp] > sbest)
			{
				_h = tmp;
				sbest = m_CTable[_h];
				continue;
			}
		}
		_c = sbest;
	}
	void getTotal(ITYPE& _v)
	{
		_v = m_Total;
	}
	void getScoresLine(const size_t& _nbTargets, std::string& _line)
	{
		std::stringstream ss;
                ITYPE num = 0;
		for(size_t t = 0; t < _nbTargets; t++)
                {
			num = m_ITable[t]==m_Token?m_CTable[t]:0;
                        ss << "," << num;
                }
                string result = ss.str();
		_line = result;
	}

	private:

	ITYPE 	m_ITable[MAXLEN];
	ITYPE	m_CTable[MAXLEN];
	ITYPE	m_Indexes[MAXLEN];
	ITYPE   m_Token;
        ITYPE   m_IBest;
	ITYPE 	m_CBest;
	ITYPE	m_Total;
	ILBL	m_ISize;
};

#endif
