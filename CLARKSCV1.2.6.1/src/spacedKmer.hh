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

#ifndef SPACEDKMER_HH
#define SPACEDKMER_HH

#include <cstdlib>
#include <iomanip>
#include <vector>
#include <cstring>
#include <stdint.h>

class spacedKmer
{
	private:
	std::string	m_name;
	std::string	m_folder;
	size_t		m_weight;
	size_t		m_length;
	bool 		m_mask[32];
	int		m_table[256];
	
	public:
	spacedKmer():  m_weight(22), m_length(31), m_name("T59575"), m_folder("T59575")
	{
		init();
		m_mask[1] = false;
		m_mask[2] = false;
                m_mask[3] = false;
		m_mask[4] = false;
		m_mask[5] = false;
                m_mask[6] = false;
                m_mask[7] = false;
		m_mask[8] = false;
		m_mask[16] = false;
		//print();
	}
	spacedKmer(const std::string& _name): m_name(_name), m_folder(_name)
	{
		string _mask;
		if (_name =="T295")     {_mask="1111*111*111**1*111**1*11*11111"; 	}
		if (_name =="T38570")   {_mask="11111*1*111**1*11*111**11*11111";       }
		if (_name =="T58570")   {_mask="11111*1**111*1*11*11**111*11111";       }

		init();
		m_weight = _mask.size();
                m_length = _mask.size();
                //m_mask.resize(m_length, true);
                for(size_t t = 0; t < m_length; t++)
                {       if (_mask[t] == '0' || _mask[t] == '*')
                        {       m_weight--;
                                m_mask[t] = false;      }
                }
                //print();
	}

	spacedKmer(const std::string& _name, const std::string& _mask): m_name(_name), m_folder(_name)
	{
		init();
		m_weight = _mask.size();
		m_length = _mask.size();
		//m_mask.resize(m_length, true);
		for(size_t t = 0; t < m_length; t++)
		{	if (_mask[t] == '0' || _mask[t] == '*')
			{	m_weight--;
				m_mask[t] = false;	}
		}
		//print();
	}
	~spacedKmer(){}
	std::string getFolder() const 
	{	return m_folder;}

	std::string getName() const
	{	return m_name;	}

	size_t getWeight() const
        {       return m_weight;}	

	size_t getLength() const
        {       return m_length;}

	bool isFwdValid(const uint8_t* _map, const size_t& _i, const size_t& _max, uint64_t& _kmf) const
	{
		if (_max - _i < m_length)
                {
		        return false;
		}
		size_t c = 0, i_c = _i;
		_kmf = 0;
		while (true)
		{
			if (c == m_length)
			{	return true;	}
			if (i_c == _max && c < m_length)
                        {	return false;	}
			if (m_table[_map[i_c]] == -10)
			{	i_c++;
				continue;
			}
			if (m_table[_map[i_c]] < 0)
                        {       
                                return false;
                        }
			if (m_table[_map[i_c]] != 4)
			{
				_kmf <<= 2;
				_kmf ^= (m_table[_map[i_c]]);
				i_c++;
				c++;
				continue;
			}
			if (!m_mask[c])
			{
				_kmf <<= 2;
                                _kmf ^= 0;
				c++;
				i_c++;	
				continue;
			}
			return false;	
		}
	}
	bool isRvsValid(const uint8_t* _map, const size_t& _i, const size_t& _max, uint64_t& _kmr) const
	{
		if (_max - _i < m_length)
                {
		        return false;
		}
                size_t c = 0, i_c = _i;
                _kmr = 0;
		while (true)
                {
                        if (c == m_length)
                        {       return true;    }       
                        if (i_c == _max && c < m_length)
                        {       return false;   }
                        if (m_table[_map[i_c]] == -10)
                        {       i_c++;
                                continue;
                        }
                        if (m_table[_map[i_c]] < 0)
                        {
                                return false;
                        }
                        if (m_table[_map[i_c]] != 4)
                        {
                                _kmr <<= 2;
                                _kmr ^= (m_table[_map[i_c]]);
                                i_c++;
                                c++;
                                continue;
                        }
                        if (!m_mask[m_length - 1 - c])
                        {       
                                _kmr <<= 2;
                                _kmr ^= 0;
                                c++;
                                i_c++; 
                                continue;
                        }
                        return false;
                }
	}
	private:
	void init()
	{
		for(size_t u = 0; u < 256; u++)
                {       
                        m_table[u] = -1;
                }
                m_table['A']  = 3; m_table['C'] = 2; m_table['G'] = 1; m_table['T'] = 0; m_table['U'] = 0;
                m_table['a']  = 3; m_table['c'] = 2; m_table['g'] = 1; m_table['t'] = 0; m_table['u'] = 0;
                m_table['\n'] = -10;
                m_table['n']  = 4; m_table['N'] = 4;

		m_table['M']  = 4; m_table['R'] = 4; m_table['W'] = 4; m_table['V'] = 4; m_table['D'] = 4;
                m_table['K']  = 4; m_table['Y'] = 4; m_table['S'] = 4; m_table['H'] = 4; m_table['B'] = 4;
                m_table['m']  = 4; m_table['r'] = 4; m_table['w'] = 4; m_table['v'] = 4; m_table['d'] = 4;
                m_table['k']  = 4; m_table['y'] = 4; m_table['s'] = 4; m_table['h'] = 4; m_table['b'] = 4;

		for(size_t t = 0; t < 32 ; t++)
		{
			m_mask[t] = true;
		}
	}
	void print() const
	{		
		cerr << "Spaced seed info: " << m_name << " " << m_weight <<  " " << m_length << endl;
		for(size_t t = 0; t < m_length; t++)
		{	cerr << (m_mask[t]?"1":"*") ;	}
		cerr << endl;
	}
};

#endif
