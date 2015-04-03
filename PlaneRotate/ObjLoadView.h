//This example program is created by thecplusplusuy for demonstration purposes. It's a simple 3D model loader (wavefront (.obj)), which is capable to load materials and UV textures:
//http://www.youtube.com/user/thecplusplusguy
//Free source, modify if you want, LGPL licence (I guess), I would be happy, if you would not delete the link
//so other people can see the tutorial
//this file is the objloader.h
#ifndef OBJLOADER_H
#define OBJLOADER_H


#include <SDL.h>
#include <gl/glew.h>
#include <gl/freeglut.h>
#include <gl/GLAUX.H>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <assert.h>
#include <memory.h>
#include  "gdiam.h"

using namespace std;
//This struct contain 3 floats and a constructor, it's used for vertexes and normal vectors
struct coordinate
{
	float x, y, z;
	coordinate(float a, float b, float c);
};

//This structure is store every property of a face
struct face
{
	int facenum;    //the number of the face (it's start from 1 not 0, so if you use it as an index, subtract 1 from it), it's used for the normal vectors
	bool four;              //if true, than it's a quad else it's a triangle
	int faces[4];   //indexes for every vertex, which makes the face (it's start from 1 not 0, so if you use it as an index, subtract 1 from it)
	int texcoord[4];        //indexes for every texture coorinate that is in the face (it's start from 1 not 0, so if you use it as an index, subtract 1 from it)
	int normal[4];
	int mat;                                        //the index for the material, which is used by the face
	face(int facen, int f1, int f2, int f3, int t1, int t2, int t3, int n1, int n2, int n3,int m);        //constuctor for triangle
	face(int facen, int f1, int f2, int f3, int f4, int t1, int t2, int t3, int t4, int n1, int n2, int n3, int n4,int m);  //-"- for quad
};

//this is a structure, which contain one material
struct material
{
	std::string name;       //the name of the material
	float alpha, ns, ni;      //some property, alpha, shininess, and some other, which we not used
	float dif[3], amb[3], spec[3];    //the color property (diffuse, ambient, specular)
	int illum;      //illum - we not use it
	int texture;    //the id for the texture, if there is no texture than -1
	material(const char* na, float al, float n, float ni2, float* d, float* a, float* s, int i, int t);
};

//texture coorinate (UV coordinate), nothing to explain here
struct texcoord
{
	float u, v;
	texcoord(float a, float b);
};

//the main class for the object loader
class objloader
{
	GBBox   bbx;
	gdiam_real  * points;
public:
	GLfloat x_min ;
	GLfloat y_min ;
	GLfloat z_min ;

	GLfloat x_max ;
	GLfloat y_max ;
	GLfloat z_max ;

	GLfloat delt_width;//grid分割小方块的宽度，分割的方块为正方形

	int divide_num;//grid边界分割基数

public:
	vector<vector<face>> faces_in_grid;
	std::vector<std::string*> coord;        //every line of code from the obj file
	std::vector<coordinate*> vertex;        //all vertexes
	std::vector<face*> faces;                                       //all faces
	std::vector<coordinate*> normals;       //all normal vectors
	std::vector<unsigned int> texture;//the id for all the textures (so we can delete the textures after use it)
	std::vector<unsigned int> lists;        //the id for all lists (so we can delete the lists after use it)
	std::vector<material*> materials;       //all materials
	std::vector<texcoord*> texturecoordinate;       //all texture coorinate (UV coordinate)
	bool ismaterial, isnormals, istexture;    //obvious
	unsigned int loadTexture(const char* filename); //private load texture function
	
public:
	objloader();
	~objloader();   //free the textures and lists
	void clean();   //free all of the used memory
	int load(const char* filename); //the main model load function
	void CalculateBoundingBox();//计算boundingbox的最大最小坐标
	void DrawBoundingBox();//calculate model bounding box---obb,aabb
	void DrawGridCell();//绘制grid cell
	void DividFace2Grid();//对每个cell中的面片进行统计
	void objloader::DrawGridCellFace(const vector<face> &grid_cell_faces_in);//绘制cell中的面片
};

#endif