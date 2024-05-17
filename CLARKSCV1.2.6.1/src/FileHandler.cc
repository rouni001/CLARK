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

#include "./FileHandler.hh"
#include "./file.hh"
using namespace std;

FileHandler::FileHandler(const char* filename,const int& _nbCPU, const size_t& _maxNbReads):
        _map(NULL),
        i_Pos(nbCPU,0),
        _filename(filename),
        nbCPU(_nbCPU),
        max(nbCPU,0),
        maxNbReads(_maxNbReads),
        _fileSize(0),
        _idxFrag(0),
	i_PosDone(nbCPU,false),
	_rIndex(_nbCPU,0)
{}

FileHandler::~FileHandler()
{
	free(_map);
	_map = NULL;
	Close();
}

bool	FileHandler::Next()
{
	if (_idxFrag >= maxPos.size() -1)
		return false;
	size_t sizeToLoad = maxPos[_idxFrag+1] - maxPos[_idxFrag];
	if  (sizeToLoad == 0)
	{
		cerr << "Problem with the memory size to allocate."<< endl;
		exit(1);
	}
	if (_idxFrag == 0)
	{
		_map = (uint8_t *) calloc(sizeToLoad, 1);
	}
	else
	{
		_map = (uint8_t *) realloc(_map, sizeToLoad);
	}

	if (_map == NULL)
	{
		cerr << "Failed to allocate the memory space needed."<< endl;
		exit(1);
	}
	size_t nbRead  = fread(_map, 1, sizeToLoad, fd);
	if (nbRead != sizeToLoad)
	{
		cerr << "Issue to load all characters (" << nbRead << " B read instead of " << sizeToLoad<< endl; 
		exit(1);
	}	
	i_PosDone.clear();
	i_PosDone.resize(nbCPU, false);
	
	i_Pos.clear();
	i_Pos.resize(nbCPU,0);

	_rIndex.clear();
	_rIndex.resize(nbCPU,0);

	SetPositions();

	_idxFrag++;
	return true;
}

bool 	FileHandler::Open()
{
	return false;
}

bool 	FileHandler::Close()
{
	if (fd != NULL)
	{
		fclose(fd);
		fd = NULL;
		return true;
	}
	return false;
}

size_t	FileHandler::GetNbCPU() const
{	return nbCPU;}

uint64_t FileHandler::GetReadsCount() const
{
	return _nbReads;
}

uint64_t FileHandler::GetCurrReadsCount() const
{
	if (_fileSize == 0 && _idxFrag == 0)
	{	
		cerr << "There is no read loaded. " << endl;
		return 0;	
	}
	return currNbReads[_idxFrag]-currNbReads[_idxFrag-1];
}

uint64_t FileHandler::GetReadID(const int& i_cpu) const
{
        return _rIndex[i_cpu];
}

bool FileHandler::GetRead(const int& i_cpu, std::string& out, std::string& id)
{
	return false;
}

bool FileHandler::GetRead(const int& i_cpu, uint8_t* out, uint32_t& size, std::string& id)
{
	return false;
}

bool 	FileHandler::SetPositions()
{
	return false;
}	

bool FileHandler::isStart() const
{
	return _idxFrag == 1;
}

bool FileHandler::isEnd() const
{
        return _idxFrag == maxPos.size() -1;
}

bool 	FileHandler::Getline(const int& i_cpu, std::string& out, const bool isStoring)
{
	out = ""; //out.clear();
        if (i_PosDone[i_cpu])
                return false;
        if (isStoring)
        {
                uint64_t i = 0;
                for(i = i_Pos[i_cpu]; i < max[i_cpu]; i++)
                {
                        if (_map[i] == 10)
                                break;
                        out.push_back((char) _map[i]);
                }
                i_Pos[i_cpu] = i+1;
                if (i_Pos[i_cpu] >= max[i_cpu])
                {
                        i_PosDone[i_cpu] = true;
                }
                return true;
        }
        uint64_t i = 0;
        for(i = i_Pos[i_cpu]; i < max[i_cpu] && _map[i] != 10; i++)
        {}
        i_Pos[i_cpu] = i+1;
        if (i_Pos[i_cpu] >= max[i_cpu])
        {
                i_PosDone[i_cpu] = true;
        }
        return true;
}

bool 	FileHandler::isOver(const int& i_cpu) const
{
	return i_PosDone[i_cpu];
}

uint64_t FileHandler::Size() const
{
	return _fileSize;
}

