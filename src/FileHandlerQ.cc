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
 * @project: CLARK, Metagenomic and Genomic Sequences Classification project.
 * @note: C++ IMPLEMENTATION supported on latest Linux and Mac OS.
 *
 */

#include <iostream>
#include <cstdlib>

#include "./FileHandlerQ.hh"
#include "./file.hh"
using namespace std;

FileHandlerQ::FileHandlerQ(const char* filename,const int& _nbCPU, const size_t& _maxNbReads): FileHandler(filename,_nbCPU, _maxNbReads)
{}

FileHandlerQ::~FileHandlerQ()
{
}

bool FileHandlerQ::GetRead(const int& i_cpu, std::string& out, std::string& id)
{
        id = "";
	out = ""; 
        if (i_PosDone[i_cpu])
                return false;
	
        bool stop = false;
	
	uint64_t i = i_Pos[i_cpu]+1;
        while (i < max[i_cpu] && _map[i] != 10)
        {
		stop = stop || _map[i] == ' ' || _map[i] =='\t';
		if (!stop)
			id.push_back(_map[i]);
		i++;
	}
	i++;
        while (i < max[i_cpu] && _map[i] != 10)
        {
                out.push_back((char) _map[i++]);
        }
        i++;
        while (i < max[i_cpu] && _map[i++] != 10)
        {}
        while (i < max[i_cpu] && _map[i++] != 10)
        {}
        i_Pos[i_cpu] = i;

        if (i_Pos[i_cpu] >= max[i_cpu])
        {
                i_PosDone[i_cpu] = true;
        }
	_rIndex[i_cpu]++;
        return true;
}

bool FileHandlerQ::GetRead(const int& i_cpu, uint8_t* out, uint32_t& size, std::string& id)
{
	size = 0;
	id = "";
        if (i_PosDone[i_cpu])
                return false;
	if (_map[i_Pos[i_cpu]] != '@')
	{
		cerr << "Failed to read from the header (" << i_cpu<< "), read number: " << _rIndex[i_cpu] << endl;
		exit(1);
	}
        bool stop = false;
        uint64_t i = i_Pos[i_cpu]+1;
        while (i < max[i_cpu] && _map[i] != 10)
        {
                stop = stop || _map[i] == ' ' || _map[i] =='\t';
                if (!stop)
                        id.push_back(_map[i]);
                i++;
        }
        i++;
        while (i < max[i_cpu] && _map[i] != 10)
        {
                out[size++] = ((uint8_t) _map[i++]);
        }
        i++;
        while (i < max[i_cpu] && _map[i++] != 10)
        {}
        while (i < max[i_cpu] && _map[i++] != 10)
        {}
        i_Pos[i_cpu] = i;

        if (i_Pos[i_cpu] >= max[i_cpu])
        {
                i_PosDone[i_cpu] = true;
        }
        _rIndex[i_cpu]++;
	return true;
}


bool 	FileHandlerQ::Open()
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
	while (true)
        {
                t = 0 ;
                for(t = 0; t < 4 && getLineFromFile(fd, line); t++)
                {       cumulSize += line.size()+1;     }
                if (t == 0)
                        break;
                if (t < 4)
                {       cerr << "Bad format file!"<< endl;
                        exit(1);
                }
                cumulNbReads++;
                _nbReads++;
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
		if (cumulSize > maxPos[maxPos.size()-1] + 100 * nbCPU)
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

bool 	FileHandlerQ::SetPositions()
{
        uint64_t fragmentSize = maxPos[_idxFrag+1]-maxPos[_idxFrag];
        uint64_t bigSteps = fragmentSize/nbCPU;
        uint64_t _Size = fragmentSize;

        posReads.resize(nbCPU, fragmentSize);
        posReads[0] = 1;
	i_Pos[0] = 0;
        if (fragmentSize < 100*nbCPU)
        {
		i_PosDone[0] = false;
                i_Pos[0] = 0;
                max[0] = _Size;
                for(size_t t = 1; t < nbCPU; t++)
                {
                        i_Pos[t] = _Size;
                        max[t] = _Size;
			i_PosDone[t] = true;
                }
		_rIndex.resize(nbCPU,0);
                return true;
        }

        int Letter[256];
        for(size_t t = 0; t < 256 ; t++)
        {       Letter[t] = -1; }

        // Letter definition (Fastq files)
        for(size_t t = 65; t < 91; t++)
        {       Letter[t] = 0;  }
        for(size_t t = 97; t < 123; t++)
        {       Letter[t] = 0;  }
        for(size_t i_r = 1; i_r < nbCPU ; i_r++)
        {
                size_t pos1 = 0, pos2 = 0, pos3 = 0, pos4 = 0, pos5 = 0, pos6 = 0;
                size_t i = bigSteps * i_r;

                while (i < _Size && _map[i++] != '\n')
                {       }
                pos1 = i;
                while (i < _Size && _map[i++] != '\n')
                {       }
                pos2 = i;
                while (i < _Size && _map[i++] != '\n')
                {       }
                pos3 = i;
                while (i < _Size && _map[i++] != '\n')
                {       }
                pos4 = i;
                while (i < _Size && _map[i++] != '\n')
                {       }
                pos5 = i;
                while (i < _Size && _map[i++] != '\n')
                {       }
                pos6 = i;
                if (_map[pos1] == '@')
                {
                        i = pos2;
                        while (Letter[_map[i++]] >= 0)
                        {}
                        if (i == pos3 && _map[i] == '+')
                        {
                                posReads[i_r] = pos1+1;
                                continue;
                        }
                }
                if (_map[pos2] == '@')
                {
                        i = pos3;
                        while (Letter[_map[i++]] >= 0)
                        {}
                        if (i == pos4 && _map[i] == '+')
                        {       posReads[i_r] = pos2+1;
                                continue;
                        }
                }
                if (_map[pos3] == '@')
                {
                        i = pos4;
                        while (Letter[_map[i++]] >= 0)
                        {}
                        if (i == pos5 && _map[i] == '+')
                        {       posReads[i_r] = pos3+1;
                                continue;
                        }
                }
                if (_map[pos4] == '@')
                {
                        i = pos5;
                        while (Letter[_map[i++]] >= 0)
                        {}
                        if (i == pos6 && _map[i] == '+')
                        {       posReads[i_r] = pos4+1;
                                continue;
                        }
                }
        }
        for(size_t t = 0; t < nbCPU; t++)
        {
                i_Pos[t] = posReads[t]-1;
                max[t] = (t<nbCPU-1)?(posReads[t+1]-1):_Size;
        }

	_rIndex.resize(nbCPU,0);
        for(size_t i_r = 1; i_r < nbCPU ; i_r++)
        {
                size_t i = i_Pos[i_r-1];
                while (i < i_Pos[i_r])
                {
                        if (_map[i++] == '@')
                        {       _rIndex[i_r]++;
                                // @header
                                while (_map[i++] != 10)
                                {}
                                // read
                                while (_map[i++] != 10)
                                {}
                                // +header
                                while (_map[i++] != 10)
                                {}
                                // quality
                                while (_map[i++] != 10)
                                {}
                        }
                }
		if (i_r < nbCPU -1)
                	_rIndex[i_r+1] = _rIndex[i_r];
        }	
	return true;
}	

