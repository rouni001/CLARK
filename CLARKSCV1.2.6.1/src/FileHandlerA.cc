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

#include <iostream>
#include <cstdlib>

#include "./FileHandlerA.hh"
#include "./file.hh"
using namespace std;

FileHandlerA::FileHandlerA(const char* filename,const int& _nbCPU, const size_t& _maxNbReads):FileHandler(filename,_nbCPU,_maxNbReads)
{}

FileHandlerA::~FileHandlerA()
{
}

bool FileHandlerA::GetRead(const int& i_cpu, std::string& out, std::string& id)
{
	id = "";
	out = ""; 
	if (i_PosDone[i_cpu])
		return false;

	uint64_t i = i_Pos[i_cpu]+1;
	bool stop = false;
	while (i < max[i_cpu] && _map[i] != 10)
	{
		stop = stop || (_map[i] == ' ' || _map[i] =='\t');
		if (!stop)
			id.push_back(_map[i]);
		i++;
	}
	while (i < max[i_cpu] && _map[i] != '>')
	{
		if (_map[i] != 10)
			out.push_back((char) _map[i]);
		i++;
	}
	i_Pos[i_cpu] = i;

	if (i_Pos[i_cpu] >= max[i_cpu])
	{
		i_PosDone[i_cpu] = true;
	}
	_rIndex[i_cpu]++;
	return true;	
}

bool FileHandlerA::GetRead(const int& i_cpu, uint8_t* out, uint32_t& size, std::string& id)
{
	size = 0;
	id = "";
	if (i_PosDone[i_cpu])
		return false;

	uint64_t i = i_Pos[i_cpu]+1;
	bool stop = false;

	while (i < max[i_cpu] && _map[i] != 10)
	{
		stop = stop || _map[i] == ' ' || _map[i] =='\t';
		if (!stop)
			id.push_back(_map[i]);
		i++;
	}
	while (i < max[i_cpu] && _map[i] != '>')
	{
		if (_map[i] != 10)
		{
			out[size++] = ((uint8_t) _map[i]);
		}
		i++;
	}
	i_Pos[i_cpu] = i;

	if (i_Pos[i_cpu] >= max[i_cpu])
	{
		i_PosDone[i_cpu] = true;
	}
	_rIndex[i_cpu]++;
	return true;
}


bool 	FileHandlerA::Open()
{
	_fileSize = 0;

	fd = fopen(_filename, "r");
	if (fd == NULL)
	{
		return false;
	}
	std::string line = "";
	uint64_t cumulSize = 0, cumulNbReads = 0;

	maxPos.push_back(0);
	currNbReads.push_back(0);

	size_t t = 0;
	_nbReads = 0;

	bool done = !getLineFromFile(fd, line);
	if (line[0] != '>')
	{
		cerr << "Failed to recognize this format file."<< endl;
		return false;
	}

	while (true)
	{
		if (done)
			break;

		cumulSize += line.size()+1;
		cumulNbReads++;
		_nbReads++;

		while (true)
		{
			done = !getLineFromFile(fd, line);
			if (done)
				break;
			if (line[0] == '>')
			{
				break;
			}
			cumulSize += line.size()+1;
		}
		if (done)
			break;
		if (cumulNbReads >= maxNbReads)
		{
			maxPos.push_back(cumulSize);
			currNbReads.push_back(_nbReads);
			cumulNbReads = 0;
		}
	}
	if (2*nbCPU >= _nbReads)
        {
                nbCPU = 1;
                max.clear(); 
                max.resize(1,0),
                i_PosDone.clear();
                i_PosDone.resize(1,false);
                _rIndex.clear();
                _rIndex.resize(1,0);
        }
	if (cumulSize > maxPos[maxPos.size()-1])
	{
		if (cumulSize >= maxPos[maxPos.size()-1] + 100 * nbCPU)
		{	
			maxPos.push_back(cumulSize);
			currNbReads.push_back(_nbReads);
		}
		else
		{
			maxPos[maxPos.size()-1] = cumulSize;
			currNbReads[maxPos.size()-1] = _nbReads;
		}
	}
	_fileSize = cumulSize;
	fclose(fd);
	fd = fopen(_filename, "r");
	return true;
}

bool 	FileHandlerA::SetPositions()
{
	uint64_t fragmentSize = maxPos[_idxFrag+1]-maxPos[_idxFrag];
	uint64_t bigSteps = fragmentSize/nbCPU;
	uint64_t _Size = fragmentSize;

	posReads.resize(nbCPU, fragmentSize);
	_rIndex.resize(nbCPU,0);

	for(size_t i_r = 1; i_r < nbCPU ; i_r++)
	{
		size_t i = bigSteps * (i_r-1);
		while (i < bigSteps * i_r)
		{
			if (_map[i++] == '>')
			{	
				_rIndex[i_r]++;	
				while ( _map[i++] != '\n' )
				{}
			}
		}
		if (i_r < nbCPU - 1)
		{	_rIndex[i_r+1] = _rIndex[i_r];	}
	}

	posReads[0] = 1;
	i_Pos[0] = 0;
	if (_map[0] != '>')
	{
		cerr << "Failed to divide properly the fasta file."<< endl;
		exit(1);
	}

	if (fragmentSize < 100*nbCPU)
	{
		i_Pos[0] = 0;
		max[0] = _Size;
		for(size_t t = 1; t < nbCPU; t++)
		{
			i_Pos[t] = _Size;
			max[t] = _Size;
			i_PosDone[t] = true;
		}
		return true;
	}
	for(size_t i_r = 1; i_r < nbCPU ; i_r++)
	{
		size_t i = bigSteps * i_r;
		while (i < _Size && _map[i++] != '>')
		{}
		posReads[i_r] = i;
	}
	for(size_t t = 0; t < nbCPU; t++)
	{
		i_Pos[t] = posReads[t]-1;
		max[t] = (t<nbCPU-1)?(posReads[t+1]-1):_Size;
	}

	return true;
}	
