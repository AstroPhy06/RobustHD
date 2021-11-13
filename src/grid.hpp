#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <iostream>


template <typename T>
class Grid
{
private:
	double _Lx;	// metric size in x
	double _Ly;	// metric size in y
	double _Nx;	// number of elements in x
	double _Ny;	// number of elements in y
	
	double _dx;
	double _dy;
	
	std::vector<T> _data;
	
public:

	Grid(double Lx, double Ly, double Nx, double Ny): _Lx(Lx), _Ly(Ly), _Nx(Nx), _Ny(Ny), _data((Nx+3)*(Ny+3))
	{
		_dx = _Lx/_Nx;
		_dy = _Ly/_Ny;
	}
	
	/// \brief get reference to value of the physical grid (without dead zones)
	inline T& at(int x, int y)
	{
		return _data[(2+y)*(_Nx+3)+2+x];
	}
	
	/// \brief get value of the physical grid (without dead zones)
	inline T get_at(int x, int y)
	{
		return _data[(2+y)*(_Nx+3)+2+x];
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////
	
	std::vector<T> row(int i)
	{
	   std::vector<T> newVec(_data.begin() + (i+2)*(_Nx+3), _data.begin() + (i+3)*(_Nx+3));
	   return newVec;
	}
	
	void setrow(int i, std::vector<T> r)
	{
	   for(int j=0; j<r.size(); j++)
 	   {
              _data[j+(i+2)*(_Nx+3)] = r[j];
           }
	}
	
	std::vector<T> col(int i)
	{
	   std::vector<T> newVec(_Ny+3);
	   for(int j=0; j<newVec.size(); j++)
 	   {
 	      newVec[j] = get_at(i,j-2);
 	   }
	   return newVec;
	}
	
	void setcol(int i, std::vector<T> r)
	{
	   for(int j=0; j<r.size(); j++)
 	   {
              _data[i+2+j*(_Nx+3)] = r[j];
           }
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////
	
	/// \brief set value of the physical grid (without dead zones)
	inline void set_at(int x, int y, T a)
	{
		_data[(2+y)*(_Nx+3)+2+x] = a;
	}
	
	void init_border_to_value(const T &val)
	{
		for(int i=0; i<_Nx+3; i++)
		{
			_data[i] = val;
			_data[i+_Nx+3] = val;
			_data[i+(_Ny+2)*(_Nx+3)] = val;
		}
		for(int j=0; j<_Ny+3; j++)
		{
			_data[j*(_Nx+3)] = val;
			_data[j*(_Nx+3)+1] = val;
			_data[j*(_Nx+3)+_Nx+2] = val;
		}
	}
	
	void print(int k)
	{
	  for (int i=0; i<_Ny+3; i++)
	  {
	    for (int j=0; j<_Nx+3; j++)
	    {
	      std::cout<<_data[j+i*(_Nx+3)][k]<<"	";
	    }
	    std::cout<<std::endl;
	  }
	  std::cout<<std::endl;
	}
	
	size_t size(){return _data.size();}
	
	int Nx(){return _Nx;}
	int Ny(){return _Ny;}
	
	double dx(){return _dx;}
	double dy(){return _dy;}
	
	void swap(Grid &g) {g._data.swap(_data);}
	void swap(std::vector<T> &v) {v.swap(_data);}
	
	std::vector<T> data(){return _data;}
	
};



#endif
