#include <Eigen/Dense>
#include "voxelize.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <cxxopts.hpp>
#include <cuda_profiler_api.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define BLOCK 512
#define MAXLENGTH 256

using namespace Eigen;
using namespace std;

bool verbose_output = false, parallel = false;
int resolution = 16;	
float voxelsize = 100.0f;

Vector3f min_corner, max_corner;

vector<tetrahedra> mesh;
vector<Eigen::Vector3f> vertex;
vector<float> sample_coord_x, sample_coord_y, sample_coord_z;
bool voxels[MAXLENGTH][MAXLENGTH][MAXLENGTH] = {0};
int MAXSIZE = MAXLENGTH * MAXLENGTH * MAXLENGTH;

__device__
int GetIndex(int i, int j, int k)
{
	return (i * MAXLENGTH + j ) * MAXLENGTH + k;
}

inline float min(float a, float b, float c, float d)
{
	return min(min(min(a,b),c),d);
}

inline float max(float a, float b, float c, float d)
{
	return max(max(max(a,b),c),d);
}

void NodesInput(const char *fname)
{
	int input_count, id, polyface = 3, type1 = 0, type2 = 0;
	float x_coord, y_coord, z_coord;

	ifstream fin(fname);
	fin >> input_count >> polyface >> type1 >> type2;

	assert(polyface == 3 && type1 == 0 && type2 == 0);
	vertex.resize(input_count);

	for(int i = 1; i <= input_count; ++i)
	{
		fin >> id >> x_coord >> y_coord >> z_coord;
		assert(id == i);
		vertex[i] << x_coord, y_coord, z_coord;
	}

	if(verbose_output)
		cout << "NodeCount: " << id << endl;
	string message;
	fin.ignore();
	getline(fin, message);
	if(verbose_output)
		cout << message << endl;
}

void ElementInput(const char *fname)
{
	int input_count, id;
	int poly_face = 4, type = 0, id_A, id_B, id_C, id_D;
	
	ifstream fin(fname);
	fin >> input_count >> poly_face >> type;
	
	assert(poly_face == 4 && type == 0);
	mesh.resize(input_count + 1);

	for(int i = 1; i <= input_count; ++i)
	{
		fin >> id >> id_A >> id_B >> id_C >> id_D;
		assert(id == i);
		tetrahedra t;
		t.a = vertex[id_A];
		t.b = vertex[id_B];
		t.c = vertex[id_C];
		t.d = vertex[id_D];
	
		t.bound_min << min(t.a(0), t.b(0), t.c(0), t.d(0)), min(t.a(1), t.b(1), t.c(1), t.d(1)), min(t.a(2), t.b(2), t.c(2), t.d(2));
		t.bound_max << max(t.a(0), t.b(0), t.c(0), t.d(0)), max(t.a(1), t.b(1), t.c(1), t.d(1)), max(t.a(2), t.b(2), t.c(2), t.d(2));
	
		mesh[i] = t;
	}

	if(verbose_output)
		cout << "ElementCount: " << id << endl;
	string message;
	fin.ignore();
	getline(fin, message);
	if(verbose_output)
		cout << message << endl;
}

void MshInput(string fname = "")
{
	string nd = "$Nodes", ednd = "$EndNodes", ele = "$Elements", edele = "$EndElements";
	string s;
	int input_count;
	int id;
	float x_coord, y_coord, z_coord;
	int poly_face = 4, type = 0, id_A, id_B, id_C, id_D;

	if(fname.empty())
	{
		while(cin >> s)
		{
			if(s == nd)
			{
				cin >> input_count;
				vertex.resize(input_count);
				for(int i = 1; i <= input_count; ++i)
				{
					cin >> id >> x_coord >> y_coord >> z_coord;
					assert(id == i);
					vertex[i] << x_coord, y_coord, z_coord;
				}

				if(verbose_output)
					cout << "NodeCount: " << id << endl;
				cin >> s;
				assert(s == ednd);
			}
			if(s == ele)
			{
				cin >> input_count;
				mesh.resize(input_count);
				for(int i = 1; i <= input_count; ++i)
				{
					cin >> id >> poly_face >> type >> id_A >> id_B >> id_C >> id_D;
					assert(id == i && poly_face == 4 && type == 0);
					tetrahedra t;
					t.a = vertex[id_A];
					t.b = vertex[id_B];
					t.c = vertex[id_C];
					t.d = vertex[id_D];
				
					t.bound_min << min(t.a(0), t.b(0), t.c(0), t.d(0)), min(t.a(1), t.b(1), t.c(1), t.d(1)), min(t.a(2), t.b(2), t.c(2), t.d(2));
					t.bound_max << max(t.a(0), t.b(0), t.c(0), t.d(0)), max(t.a(1), t.b(1), t.c(1), t.d(1)), max(t.a(2), t.b(2), t.c(2), t.d(2));
				
					mesh[i] = t;
				}
				if(verbose_output)
					cout << "ElementCount: " << id << endl;
				cin >> s;
				assert(s == edele);
			}
		}
	}
	else
	{
		ifstream fin(fname.c_str());
		while(fin >> s)
		{
			if(s == nd)
			{
				fin >> input_count;
				vertex.resize(input_count + 1);
				for(int i = 1; i <= input_count; ++i)
				{
					fin >> id >> x_coord >> y_coord >> z_coord;
					assert(id == i);
					vertex[i] << x_coord, y_coord, z_coord;
				}

				if(verbose_output)
					cout << "NodeCount: " << id << endl;
				fin >> s;
				assert(s == ednd);
			}
			if(s == ele)
			{
				fin >> input_count;
				mesh.resize(input_count);
				for(int i = 1; i <= input_count; ++i)
				{
					fin >> id >> poly_face >> type >> id_A >> id_B >> id_C >> id_D;
					assert(id == i && poly_face == 4 && type == 0);
					tetrahedra t;
					t.a = vertex[id_A];
					t.b = vertex[id_B];
					t.c = vertex[id_C];
					t.d = vertex[id_D];
				
					t.bound_min << min(t.a(0), t.b(0), t.c(0), t.d(0)), min(t.a(1), t.b(1), t.c(1), t.d(1)), min(t.a(2), t.b(2), t.c(2), t.d(2));
					t.bound_max << max(t.a(0), t.b(0), t.c(0), t.d(0)), max(t.a(1), t.b(1), t.c(1), t.d(1)), max(t.a(2), t.b(2), t.c(2), t.d(2));
				
					mesh[i] = t;
				}
				if(verbose_output)
					cout << "ElementCount: " << id << endl;
				fin >> s;
				assert(s == edele);
			}
		}
	}
	return;
}

void NodesElementsInput(string fname)
{
	NodesInput((fname + ".node").c_str());
	ElementInput((fname + ".ele").c_str());
}

void StandardInput()
{
	int T;
	cin >> resolution;
	cin >> T;
	while(T--)
	{
		tetrahedra t;
		float inx, iny, inz;
		cin >> inx >> iny >> inz;
		t.a << inx, iny, inz;
		cin >> inx >> iny >> inz;
		t.b << inx, iny, inz;
		cin >> inx >> iny >> inz;
		t.c << inx, iny, inz;
		cin >> inx >> iny >> inz;
		t.d << inx, iny, inz;

		t.bound_min << min(t.a(0), t.b(0), t.c(0), t.d(0)), min(t.a(1), t.b(1), t.c(1), t.d(1)), min(t.a(2), t.b(2), t.c(2), t.d(2));
		t.bound_max << max(t.a(0), t.b(0), t.c(0), t.d(0)), max(t.a(1), t.b(1), t.c(1), t.d(1)), max(t.a(2), t.b(2), t.c(2), t.d(2));
		
		mesh.push_back(t);
	}
}

void BoundingBox()
{
	size_t mesh_size = mesh.size();
	#ifdef _OPENMP
	#pragma omp parallel for num_threads(8) if(parallel)
	#endif
	for(int i = 0; i < mesh_size; ++i)
	{
		min_corner << min(mesh[i].bound_min(0), min_corner(0)), min(mesh[i].bound_min(1), min_corner(1)), min(mesh[i].bound_min(2), min_corner(2));
		max_corner << max(mesh[i].bound_max(0), max_corner(0)), max(mesh[i].bound_max(1), max_corner(1)), max(mesh[i].bound_max(2), max_corner(2));
	}
	voxelsize = (max_corner - min_corner).maxCoeff() / resolution;
}

__device__
float deter(Matrix3f M)
{
	return M(0,0) * (M(1,1)*M(2,2) - M(1,2)*M(2,1)) - M(0,1) * (M(1,0)*M(2,2) - M(1,2)*M(2,0)) + M(0,2) * (M(1,0)*M(2,1) - M(1,1)*M(2,0));
}

__device__
bool Philipp(Vector3f A, Vector3f B, Vector3f C, Vector3f D, Vector3f P)
{
	Vector3f a = A - P, b = B - P, c = C - P, d = D - P;
	Matrix3f Da, Db, Dc, Dd;
	Da << b, c, d;
	Db << a, c, d;
	Dc << a, b, d;
	Dd << a, b, c;

	if((deter(Da) * deter(Dc) >= 0) && (deter(Db) * deter(Dd) >= 0) && (deter(Da) * deter(Db) <= 0))
		return true;
	return false;
}

__device__
int LocateVoxel(int reso, int tid)
{
	int yz = tid / reso;
	int z = yz / reso;
	int y = yz % reso;
	int x = tid % reso;
	return (z * MAXLENGTH + y) * MAXLENGTH + x;
}

__global__
void MeshIterKernel(int mesh_size, int resolution, bool *vox_d, tetrahedra* mesh_d, float *sample_coord_x, float *sample_coord_y, float *sample_coord_z, float voxelsize)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(tid >= mesh_size)return;

	Vector3f min_corner;
	min_corner << sample_coord_x[0], sample_coord_y[0], sample_coord_z[0];
	Vector3f start_sample = mesh_d[tid].bound_min - min_corner;
	Vector3f end_sample = mesh_d[tid].bound_max - min_corner;
	Vector3f st = start_sample / voxelsize;
	Vector3f ed = end_sample / voxelsize;

	if(tid == 0)
	{
		printf("Mesh_size: %d, Resolutoin: %d, Voxelsize: %f", mesh_size, resolution, voxelsize);
	}
	for(int i = int(st(0)); i <= ed(0); ++i)
	{
		for(int j = int(st(1)); j <= ed(1); ++j)
		{
			for(int k = int(st(2)); k <= ed(2); ++k)
			{
				if(vox_d[GetIndex(i,j,k)]) continue;
				Vector3f p;
				p << sample_coord_x[i], sample_coord_y[j], sample_coord_z[k];
				if( Philipp(mesh_d[tid].a, mesh_d[tid].b, mesh_d[tid].c, mesh_d[tid].d, p) )
				{
					vox_d[GetIndex(i,j,k)] = 1;
				}
			}
		}
	}
	return;
}

void voxelize()
{
	sample_coord_x.resize(resolution);
	sample_coord_y.resize(resolution);
	sample_coord_z.resize(resolution);
	sample_coord_x[0] = min_corner(0);
	sample_coord_y[0] = min_corner(1);
	sample_coord_z[0] = min_corner(2);
	// lookup table for the coords
	for(int i = 1; i < resolution; ++i)
	{
		sample_coord_x[i] = sample_coord_x[i - 1] + voxelsize;
		sample_coord_y[i] = sample_coord_y[i - 1] + voxelsize;
		sample_coord_z[i] = sample_coord_z[i - 1] + voxelsize;
	}

	bool *vox_d;
	tetrahedra* mesh_d;
	float *sample_coord_x_d, *sample_coord_y_d, *sample_coord_z_d;
	size_t mesh_size = mesh.size();
	size_t host_data_size = mesh_size * sizeof(tetrahedra);
	size_t sample_coord_size = resolution * sizeof(float);

	cudaMalloc((void**)&vox_d, MAXSIZE * sizeof(bool));
	cudaMalloc((void**)&mesh_d, host_data_size);
	cudaMalloc((void**)&sample_coord_x_d, sample_coord_size);
	cudaMalloc((void**)&sample_coord_y_d, sample_coord_size);
	cudaMalloc((void**)&sample_coord_z_d, sample_coord_size);

	cudaMemcpy(mesh_d, &mesh[0], host_data_size, cudaMemcpyHostToDevice);
	cudaMemcpy(sample_coord_x_d, &sample_coord_x[0], sample_coord_size, cudaMemcpyHostToDevice);
	cudaMemcpy(sample_coord_y_d, &sample_coord_y[0], sample_coord_size, cudaMemcpyHostToDevice);
	cudaMemcpy(sample_coord_z_d, &sample_coord_z[0], sample_coord_size, cudaMemcpyHostToDevice);

	int grid = (mesh_size / BLOCK) + 1;

	MeshIterKernel<<<grid,BLOCK>>>(mesh_size, resolution, vox_d, mesh_d, sample_coord_x_d, sample_coord_y_d, sample_coord_z_d, voxelsize);
	cudaDeviceSynchronize();

	cudaMemcpy(voxels, vox_d, MAXSIZE * sizeof(bool), cudaMemcpyDeviceToHost);

	cudaFree(vox_d);
	cudaFree(mesh_d);
	cudaFree(sample_coord_x_d);
	cudaFree(sample_coord_y_d);
	cudaFree(sample_coord_z_d);
}

void VoxelOutput(string fname = "")
{
	if(fname.empty())
	{
		printf("%d x %d x %d\n", resolution, resolution, resolution);
		for(int i = 0; i < resolution; ++i)
		{
			for(int j = 0; j < resolution; ++j)
			{
				for(int k = 0; k < resolution; ++k)
				{
					cout << voxels[i][j][k] << " ";
				}
			}
		}
		cout << endl;
	}
	else
	{
		ofstream fout(fname.c_str());
		fout << "Resolution: " << resolution << "x" << resolution << "x" << resolution << endl;
		for(int i = 0; i < resolution; ++i)
		{
			for(int j = 0; j < resolution; ++j)
			{
				for(int k = 0; k < resolution; ++k)
				{
					fout << voxels[i][j][k] << " ";
				}
			}
		}
		fout << endl;
	}
}

int main(int argc, char** argv)
{
	cxxopts::Options options("exe.out", "One line description of MyProgram");
	options.add_options()
	("p,parallel", "Enable parallel computation on CUDA") // a bool parameter
	("n,nodes", "Enable nodes and elements input")
  ("r,resolution", "adjust the resolution of voxel output", cxxopts::value<int>())
  ("i,input", "Input file path", cxxopts::value<std::string>())
  ("o,output", "Output file path", cxxopts::value<std::string>())
  ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
  ("h,help", "List all options")
  ;

  auto result = options.parse(argc, argv);
  
  if(result.count("help"))
  {
    cout << options.help({"", "Group"}) << endl;
    return 0;
  }

  if(result.count("parallel"))
  {
  	parallel = true;
  }

 	if(result.count("resolution"))
 	{
 		resolution = result["resolution"].as<int>();
 		cout << "Output resolution is set to " << resolution << "x" << resolution << "x" << resolution << endl;
 	}

 	if(result.count("verbose"))
 	{
 		verbose_output = result["verbose"].as<bool>();
 	}

	if(result.count("nodes"))
	{
		if(result.count("input"))
		{
			NodesElementsInput(result["input"].as<string>());
		}
		else
		{
			cout << "Nodes input requires specify input file path!\n";
			exit(0);
		}
	} 	
	else if(result.count("input"))
	{
		MshInput(result["input"].as<string>());
	}
	else
	{
		MshInput();

	}
	min_corner << INF, INF, INF;
	max_corner << -INF, -INF, -INF;
	BoundingBox();

	if(verbose_output)
		cout << min_corner << endl << endl << max_corner << endl;

  voxelize();

  if(result.count("output"))
  {
  	VoxelOutput(result["output"].as<string>());
  }
  else
  {
		VoxelOutput();
  }
	cudaError_t err = cudaGetLastError();  // add
	if (err != cudaSuccess) std::cout << "CUDA error: " << cudaGetErrorString(err) << std::endl; // add
	cudaProfilerStop();

	return 0;
}
