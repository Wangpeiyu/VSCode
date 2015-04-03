//This example program is created by thecplusplusuy for demonstration purposes. It's a simple 3D model loader (wavefront (.obj)), which is capable to load materials and UV textures:
//http://www.youtube.com/user/thecplusplusguy
//Free source, modify if you want, LGPL licence (I guess), I would be happy, if you would not delete the link
//so other people can see the tutorial
//this file is the objloader.cpp
#include "stdafx.h"
#include "ObjLoadView.h"
//nothing to explain here
coordinate::coordinate(float a, float b, float c)
{
	x = a;
	y = b;
	z = c;
}
//nothing to explain here
face::face(int facen, int f1, int f2, int f3, int t1, int t2, int t3, int n1, int n2, int n3, int m){
	facenum = facen;

	faces[0] = f1;
	faces[1] = f2;
	faces[2] = f3;

	texcoord[0] = t1;
	texcoord[1] = t2;
	texcoord[2] = t3;

	normal[0] = n1;
	normal[1] = n2;
	normal[2] = n3;
	normal[3] = 0;

	mat = m;
	four = false;
}
//nothing to explain here
face::face(int facen, int f1, int f2, int f3, int f4, int t1, int t2, int t3, int t4, int n1, int n2, int n3, int n4, int m){
	facenum = facen;
	faces[0] = f1;
	faces[1] = f2;
	faces[2] = f3;
	faces[3] = f4;

	texcoord[0] = t1;
	texcoord[1] = t2;
	texcoord[2] = t3;
	texcoord[3] = t4;

	normal[0] = n1;
	normal[1] = n2;
	normal[2] = n3;
	normal[3] = n4;

	mat = m;
	four = true;
}

//nothing to explain here
material::material(const char* na, float al, float n, float ni2, float* d, float* a, float* s, int i, int t)
{
	name = na;
	alpha = al;
	ni = ni2;
	ns = n;
	dif[0] = d[0];
	dif[1] = d[1];
	dif[2] = d[2];

	amb[0] = a[0];
	amb[1] = a[1];
	amb[2] = a[2];

	spec[0] = s[0];
	spec[1] = s[1];
	spec[2] = s[2];

	illum = i;
	texture = t;
}

//nothing to explain here
texcoord::texcoord(float a, float b)
{
	u = a;
	v = b;
}

objloader::objloader()
{
	//at default we set all booleans to false, so we don't use anything
	ismaterial = false;
	isnormals = false;
	istexture = false;
}

objloader::~objloader()
{
	//delete lists and textures
	for (std::vector<unsigned int>::const_iterator it = texture.begin(); it != texture.end(); it++)
	{
		glDeleteTextures(1, &(*it));
	}
	for (std::vector<unsigned int>::const_iterator it = lists.begin(); it != lists.end(); it++)
	{
		glDeleteLists(*it, 1);
	}
}

int objloader::load(const char* filename)
{
	std::ifstream in(filename);     //open the model file
	if (!in.is_open())
	{
		std::cout << "Nor oepened" << std::endl; //if it's not opened then error message, and return with -1
		return -1;
	}
	char buf[256];  //temp buffer
	int curmat = 0;   //the current (default) material is 0, it's used, when we read the faces
	while (!in.eof())
	{
		in.getline(buf, 256);    //while we are not in the end of the file, read everything as a string to the coord vector
		coord.push_back(new std::string(buf));
	}
	for (int i = 0; i < coord.size(); i++) //and then go through all line and decide what kind of line it is
	{
		if ((*coord[i])[0] == '#') //if it's a comment
			continue;       //we don't have to do anything with it
		else if ((*coord[i])[0] == 'v' && (*coord[i])[1] == ' ')     //if a vertex
		{
			float tmpx, tmpy, tmpz;
			sscanf(coord[i]->c_str(), "v %f %f %f", &tmpx, &tmpy, &tmpz);       //read the 3 floats, which makes up the vertex
			vertex.push_back(new coordinate(tmpx, tmpy, tmpz));       //and put it in the vertex vector
		}
		else if ((*coord[i])[0] == 'v' && (*coord[i])[1] == 'n')    //if it's a normal vector
		{
			float tmpx, tmpy, tmpz;
			sscanf(coord[i]->c_str(), "vn %f %f %f", &tmpx, &tmpy, &tmpz);
			normals.push_back(new coordinate(tmpx, tmpy, tmpz));      //basically do the same
			isnormals = true;
		}
		else if ((*coord[i])[0] == 'f')   //if it's a face
		{
			int i_face[4]; //face coorinates
			int i_t[4];    //texture coorinates
			int i_n[4];    //normal coorinates

			int isquad = count(coord[i]->begin(), coord[i]->end(), ' ');
			if ( isquad > 4)     //if this is a quad
			{
				if (coord[i]->find("//") != std::string::npos)     //if it's contain a normal vector, but not contain texture coorinate
				{
					sscanf(coord[i]->c_str(), "f %d//%d %d//%d %d//%d %d//%d", &i_face[0], &i_n[0], &i_face[1], &i_n[1], &i_face[2], &i_n[2], &i_face[3], &i_n[3]);      //read in this form
					faces.push_back(new face(i_n[3], i_face[0], i_face[1], i_face[2], i_face[3], 0, 0, 0, 0, i_n[0], i_n[1], i_n[2], i_n[3], curmat));    //and put to the faces, we don't care about the texture coorinate in this case
					//and if there is no material, it doesn't matter, what is curmat
				}
				else if (coord[i]->find("/") != std::string::npos)        //if we have texture coorinate and normal vectors
				{
					//read in this form, and put to the end of the vector
					sscanf(coord[i]->c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", &i_face[0], &i_t[0], &i_n[0], &i_face[1], &i_t[1], &i_n[1], &i_face[2], &i_t[2], &i_n[2], &i_face[3], &i_t[3], &i_n[3]);
					faces.push_back(new face(i_n[3], i_face[0], i_face[1], i_face[2], i_face[3], i_t[0], i_t[1], i_t[2], i_t[3], i_n[0], i_n[1], i_n[2], i_n[3], curmat));
				}
				else
				{
					//else we don't have normal vectors nor texture coorinate
					sscanf(coord[i]->c_str(), "f %d %d %d %d", &i_face[0], &i_face[1],& i_face[2], &i_face[3]);
					faces.push_back(new face(-1, i_face[0], i_face[1], i_face[2], i_face[3], 0, 0, 0, 0,0,0,0,0, curmat));
				}
			}//if this is a quad
			else
			{  //if it's a triangle
				//do the same, except we use one less vertex/texture coorinate/face number
				if (coord[i]->find("//") != std::string::npos)
				{
					sscanf(coord[i]->c_str(), "f %d//%d %d//%d %d//%d", &i_face[0], &i_n[0], &i_face[1], &i_n[1], &i_face[2], &i_n[2]);
					faces.push_back(new face(i_n[2], i_face[0], i_face[1], i_face[2], 0, 0, 0, i_n[0], i_n[1], i_n[2], curmat));    //and put to the faces, we don't care about the texture coorinate in this case
				}
				else if (coord[i]->find("/") != std::string::npos)
				{
					int t[3];
					sscanf(coord[i]->c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d", &i_face[0], &i_t[0], &i_n[0], &i_face[1], &i_t[1], &i_n[1], &i_face[2], &i_t[2], &i_n[2]);
					faces.push_back(new face(i_n[2], i_face[0], i_face[1], i_face[2], i_t[0], i_t[1], i_t[2], i_n[0], i_n[1], i_n[2], curmat));
				}
				else
				{
					sscanf(coord[i]->c_str(), "f %d %d %d", &i_face[0], &i_face[1], &i_face[2]);
					faces.push_back(new face(-1, i_face[0], i_face[1], i_face[2], 0, 0, 0, 0, 0, 0, curmat));
				}
			}
		}//if it's a face
		else if ((*coord[i])[0] == 'u' && (*coord[i])[1] == 's' && (*coord[i])[2] == 'e')     //use material material_name
		{
			char tmp[200];
			sscanf(coord[i]->c_str(), "usemtl %s", tmp);      //read the name of the material to tmp
			for (int i = 0; i < materials.size(); i++)     //go through all of the materials
			{
				if (strcmp(materials[i]->name.c_str(), tmp) == 0)   //and compare the tmp with the name of the material
				{
					curmat = i;       //if it's equal then set the current material to that
					break;
				}
			}
		}
		else if ((*coord[i])[0] == 'm' && (*coord[i])[1] == 't' && (*coord[i])[2] == 'l' && (*coord[i])[3] == 'l')      //material library, a file, which contain
			//all of the materials
		{
			char filen[200];
			sscanf(coord[i]->c_str(), "mtllib %s", filen);    //read the filename
			std::ifstream mtlin(filen);     //open the file
			if (!mtlin.is_open())    //if not opened error message, clean all memory, return with -1
			{
				std::cout << "connot open the material file" << std::endl;
				clean();
				return -1;
			}
			ismaterial = true;        //we use materials
			std::vector<std::string> tmp;//contain all of the line of the file
			char c[200];
			while (!mtlin.eof())
			{
				mtlin.getline(c, 200);   //read all lines to tmp
				tmp.push_back(c);
			}
			char name[200]; //name of the material
			char filename[200];     //filename of the texture
			float amb[3], dif[3], spec[3], alpha, ns, ni;        //colors, shininess, and something else
			int illum;
			unsigned int texture;
			bool ismat = false;       //do we already have a material read in to these variables?
			strcpy(filename, "\0");  //set filename to nullbyte character

			for (int i = 0; i < tmp.size(); i++) //go through all lines of the mtllib file
			{
				if (tmp[i][0] == '#')      //we don't care about comments
					continue;
				if (tmp[i][0] == 'n' && tmp[i][1] == 'e' && tmp[i][2] == 'w')  //new material
				{
					if (ismat)       //if we have a material
					{
						if (strcmp(filename, "\0") != 0)    //if we have a texture
						{
							materials.push_back(new material(name, alpha, ns, ni, dif, amb, spec, illum, texture)); //push back
							strcpy(filename, "\0");
						}
						else{
							materials.push_back(new material(name, alpha, ns, ni, dif, amb, spec, illum, -1));              //push back, but use -1 to texture             
						}
					}
					ismat = false;    //we start from a fresh material
					sscanf(tmp[i].c_str(), "newmtl %s", name); //read in the name
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'N' && tmp[i][2] == 's')      //the shininess
				{
					sscanf(tmp[i].c_str(), "%*s%f", &ns);
					ismat = true;
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'K' && tmp[i][2] == 'a')      //the ambient
				{
					sscanf(tmp[i].c_str(), "%*s %f %f %f", &amb[0], &amb[1], &amb[2]);
					ismat = true;
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'K' && tmp[i][2] == 'd')      //the diffuse
				{
					sscanf(tmp[i].c_str(), "%*s %f %f %f", &dif[0], &dif[1], &dif[2]);
					ismat = true;
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'K' && tmp[i][2] == 's')      //the specular
				{
					sscanf(tmp[i].c_str(), "%*s %f %f %f", &spec[0], &spec[1], &spec[2]);
					ismat = true;
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'N' && tmp[i][2] == 'i')      //the I don't know what is this
				{
					sscanf(tmp[i].c_str(), "%*s %f", &ni);
					ismat = true;
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'd' && tmp[i][2] == ' ')      //the alpha
				{
					sscanf(tmp[i].c_str(), "%*s %f", &alpha);
					ismat = true;
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'i' && tmp[i][2] == 'l')      //the illum (don't ask)
				{
					sscanf(tmp[i].c_str(), "%*s %d", &illum);
					ismat = true;
				}
				else if (tmp[i].length() >1 && tmp[i][1] == 'm' && tmp[i][2] == 'a')      //and the texture
				{
					if (tmp[i].find("map_Kd") != std::string::npos)
					{
						sscanf(tmp[i].c_str(), "%*s%s", filename);
						texture = loadTexture(filename);  //read the filename, and use the loadTexture function to load it, and get the id.
						ismat = true;
					}
				}
			}
			if (ismat)       //there is no newmat after the last newmat, so we have to put the last material 'manually'
			{
				if (strcmp(filename, "\0") != 0)
				{
					materials.push_back(new material(name, alpha, ns, ni, dif, amb, spec, illum, texture));
				}
				else{
					materials.push_back(new material(name, alpha, ns, ni, dif, amb, spec, illum, -1));
				}
			}
		}
		else if ((*coord[i])[0] == 'v' && (*coord[i])[1] == 't')    //back to the obj file, texture coorinate
		{
			float u, v;
			sscanf(coord[i]->c_str(), "vt %f %f", &u, &v);     //read the uv coordinate
			texturecoordinate.push_back(new texcoord(u, 1 - v));       //I push back 1-v instead of normal v, because obj file use the upper left corner as 0,0 coorinate
			//but OpenGL use bottom left corner as 0,0, so I convert it
			istexture = true;
		}
	}//get model data

	if (materials.size() == 0) //if some reason the material file doesn't contain any material, we don't have material
	{
		ismaterial = false;
	}
	else    //else we have
	{
		ismaterial = true;
	}

	std::cout << vertex.size() << " " << normals.size() << " " << faces.size() << " " << materials.size() << std::endl;     //test purposes
	//draw
	int num;
	num = glGenLists(1);      //I generate a unique identifier for the list
	glNewList(num, GL_COMPILE);
	int last = -1;    //the last material (default -1, which doesn't exist, so we use the first material)
	for (int i = 0; i < faces.size(); i++) //go throught all faces
	{
		if (last != faces[i]->mat && ismaterial)   //if we have a meterial AND the last material is not the same
		{
			//set all of the material property
			float diffuse[] = { materials[faces[i]->mat]->dif[0], materials[faces[i]->mat]->dif[1], materials[faces[i]->mat]->dif[2], 1.0 };
			float ambient[] = { materials[faces[i]->mat]->amb[0], materials[faces[i]->mat]->amb[1], materials[faces[i]->mat]->amb[2], 1.0 };
			float specular[] = { materials[faces[i]->mat]->spec[0], materials[faces[i]->mat]->spec[1], materials[faces[i]->mat]->spec[2], 1.0 };
			GLfloat matEmission[] = { 0.3, 0.1, 0.1, 1.0 };

			glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
			glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
			glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
			glMaterialf(GL_FRONT, GL_SHININESS, materials[faces[i]->mat]->ns);
			

			last = faces[i]->mat;     //set the current to last
			if (materials[faces[i]->mat]->texture == -1)       //if we don't have texture, disable it, else enable it
			{
				glDisable(GL_TEXTURE_2D);
			}
			else
			{
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, materials[faces[i]->mat]->texture); //and use it
			}
		}

		if (faces[i]->four)      //if quad
		{
			glBegin(GL_QUADS);
			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[faces[i]->facenum - 1]->x, normals[faces[i]->facenum - 1]->y, normals[faces[i]->facenum - 1]->z);    //use them
			}

			if (istexture && materials.size() > 0 && materials[faces[i]->mat]->texture != -1)  //if there are textures
			{
				//set the texture coorinate
				glTexCoord2f(texturecoordinate[faces[i]->texcoord[0] - 1]->u, texturecoordinate[faces[i]->texcoord[0] - 1]->v);
			}
			

			glVertex3f(vertex[faces[i]->faces[0] - 1]->x, vertex[faces[i]->faces[0] - 1]->y, vertex[faces[i]->faces[0] - 1]->z);

			if (istexture && materials.size() > 0 && materials[faces[i]->mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[faces[i]->texcoord[1] - 1]->u, texturecoordinate[faces[i]->texcoord[1] - 1]->v);
			}
			
			glVertex3f(vertex[faces[i]->faces[1] - 1]->x, vertex[faces[i]->faces[1] - 1]->y, vertex[faces[i]->faces[1] - 1]->z);

			if (istexture && materials.size() > 0 && materials[faces[i]->mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[faces[i]->texcoord[2] - 1]->u, texturecoordinate[faces[i]->texcoord[2] - 1]->v);
			}
			

			glVertex3f(vertex[faces[i]->faces[2] - 1]->x, vertex[faces[i]->faces[2] - 1]->y, vertex[faces[i]->faces[2] - 1]->z);

			if (istexture && materials.size() > 0 && materials[faces[i]->mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[faces[i]->texcoord[3] - 1]->u, texturecoordinate[faces[i]->texcoord[3] - 1]->v);
			}
			

			glVertex3f(vertex[faces[i]->faces[3] - 1]->x, vertex[faces[i]->faces[3] - 1]->y, vertex[faces[i]->faces[3] - 1]->z);
			glEnd();
		}
		else
		{
			glBegin(GL_TRIANGLES);
			if (isnormals)   //if there are normals
				glNormal3f(normals[faces[i]->facenum - 1]->x, normals[faces[i]->facenum - 1]->y, normals[faces[i]->facenum - 1]->z);

			if (istexture && materials.size() > 0 && materials[faces[i]->mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[faces[i]->texcoord[0] - 1]->u, texturecoordinate[faces[i]->texcoord[0] - 1]->v);
			}
			glVertex3f(vertex[faces[i]->faces[0] - 1]->x, vertex[faces[i]->faces[0] - 1]->y, vertex[faces[i]->faces[0] - 1]->z);

			if (istexture && materials.size() > 0 && materials[faces[i]->mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[faces[i]->texcoord[1] - 1]->u, texturecoordinate[faces[i]->texcoord[1] - 1]->v);
			}
			glVertex3f(vertex[faces[i]->faces[1] - 1]->x, vertex[faces[i]->faces[1] - 1]->y, vertex[faces[i]->faces[1] - 1]->z);


			if (istexture && materials.size() > 0 && materials[faces[i]->mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[faces[i]->texcoord[2] - 1]->u, texturecoordinate[faces[i]->texcoord[2] - 1]->v);
			}
			glVertex3f(vertex[faces[i]->faces[2] - 1]->x, vertex[faces[i]->faces[2] - 1]->y, vertex[faces[i]->faces[2] - 1]->z);
			
			glEnd();

		}
		
	}
	glEndList();
	lists.push_back(num);

	return num;
}

void objloader::CalculateBoundingBox()
{
	int num = vertex.size();

	points = (gdiam_point)malloc(sizeof(gdiam_point_t)* num);
	assert(points != NULL);

	for (int i = 0; i < num; i++)
	{
		points[  i * 3  ] = vertex[i]->x;
		points[i * 3 + 1] = vertex[i]->y;
		points[i * 3 + 2] = vertex[i]->z;
	}

	bbx.init();
	for (int ind = 0; ind < num; ind++)
	{
		bbx.bound(points + (ind * 3));
	}
	bbx.dump();

	free(points);
	points = NULL;

	x_min = bbx.min_coords[0];
	y_min = bbx.min_coords[1];
	z_min = bbx.min_coords[2];

	x_max = bbx.max_coords[0];
	y_max = bbx.max_coords[1];
	z_max = bbx.max_coords[2];

	divide_num = 5;

	if ((z_max - z_min)> (x_max-x_min) )//判断包围盒的长宽是否相等，因此采集样本时保证长宽相等
	{
		delt_width = (z_max - z_min) / divide_num;
	}
	else
	{
		delt_width = (x_max - x_min) / divide_num;
	}
	
	DividFace2Grid();
	
}

void objloader::DrawBoundingBox()
{
	glBegin(GL_LINES);
	glColor3f(1.0f, 1.0f, 1.0f);

	glVertex3f(x_max, y_min, z_min);
	glVertex3f(x_max, y_max, z_min);//1-2

	glVertex3f(x_max, y_max, z_min);
	glVertex3f(x_max, y_max, z_max);//2-3

	glVertex3f(x_max, y_max, z_max);
	glVertex3f(x_max, y_min, z_max);//3-4

	glVertex3f(x_max, y_min, z_max);
	glVertex3f(x_max, y_min, z_min);//4-1

	glVertex3f(x_max, y_min, z_min);
	glVertex3f(x_min, y_min, z_min);//1-8

	glVertex3f(x_max, y_max, z_min);
	glVertex3f(x_min, y_max, z_min);//2-7

	glVertex3f(x_max, y_max, z_max);
	glVertex3f(x_min, y_max, z_max);//3-6

	glVertex3f(x_max, y_min, z_max);
	glVertex3f(x_min, y_min, z_max);//4-5

	glVertex3f(x_min, y_min, z_min);
	glVertex3f(x_min, y_max, z_min);//8-7

	glVertex3f(x_min, y_min, z_min);
	glVertex3f(x_min, y_min, z_max);//8-5

	glVertex3f(x_min, y_min, z_max);
	glVertex3f(x_min, y_max, z_max);//5-6

	glVertex3f(x_min, y_max, z_max);
	glVertex3f(x_min, y_max, z_min);//6-7

	glEnd();
}

void objloader::DrawGridCell()
{
	for (int i = 1; i < divide_num; i++)//绘制横向隔断
	{
		GLfloat max_z = z_min + delt_width * i;
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);

		glVertex3f(x_min, y_max, max_z);
		glVertex3f(x_max, y_max, max_z);//上

		glVertex3f(x_min, y_min, max_z);
		glVertex3f(x_max, y_min, max_z);//下

		glVertex3f(x_min, y_max, max_z);
		glVertex3f(x_min, y_min, max_z);//左

		glVertex3f(x_max, y_max, max_z);
		glVertex3f(x_max, y_min, max_z);//右

		glEnd();
	}
	

	for (int i = 1; i < divide_num; i++)//绘制竖向隔断
	{
		GLfloat max_x = x_min + delt_width * i;
		glBegin(GL_LINES);
		glColor3f(1.0f, 1.0f, 1.0f);

		glVertex3f(max_x, y_max, z_min);
		glVertex3f(max_x, y_max, z_max);//上

		glVertex3f(max_x, y_min, z_min);
		glVertex3f(max_x, y_min, z_max);//下

		glVertex3f(max_x, y_max, z_min);
		glVertex3f(max_x, y_min, z_min);//左

		glVertex3f(max_x, y_max, z_max);
		glVertex3f(max_x, y_min, z_max);//右

		glEnd();
	}
	
}

void objloader::DrawGridCellFace(const vector<face> &grid_cell_faces_in)
{
	//draw
	int last = -1; 
	for (int i = 0; i < grid_cell_faces_in.size(); i++) //go throught all faces
	{
		if (last != (grid_cell_faces_in[i].mat) && ismaterial)   //if we have a meterial AND the last material is not the same
		{
			//set all of the material property
			float diffuse[] = { materials[grid_cell_faces_in[i].mat]->dif[0], materials[grid_cell_faces_in[i].mat]->dif[1], materials[grid_cell_faces_in[i].mat]->dif[2], 1.0 };
			float ambient[] = { materials[grid_cell_faces_in[i].mat]->amb[0], materials[grid_cell_faces_in[i].mat]->amb[1], materials[grid_cell_faces_in[i].mat]->amb[2], 1.0 };
			float specular[] = { materials[grid_cell_faces_in[i].mat]->spec[0], materials[grid_cell_faces_in[i].mat]->spec[1], materials[grid_cell_faces_in[i].mat]->spec[2], 1.0 };
			GLfloat matEmission[] = { 0.3, 0.1, 0.1, 1.0 };

			glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
			glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
			glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
			glMaterialf(GL_FRONT, GL_SHININESS, materials[grid_cell_faces_in[i].mat]->ns);


			last = grid_cell_faces_in[i].mat;     //set the current to last
			if (materials[grid_cell_faces_in[i].mat]->texture == -1)       //if we don't have texture, disable it, else enable it
			{
				glDisable(GL_TEXTURE_2D);
			}
			else
			{
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, materials[grid_cell_faces_in[i].mat]->texture); //and use it
			}
		}

		if (grid_cell_faces_in[i].four)      //if quad
		{
			glBegin(GL_QUADS);

			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[grid_cell_faces_in[i].normal[0] - 1]->x, normals[grid_cell_faces_in[i].normal[0] - 1]->y, normals[grid_cell_faces_in[i].normal[0] - 1]->z);    //use them
			}
			if (istexture && materials.size() > 0 && materials[grid_cell_faces_in[i].mat]->texture != -1)  //if there are textures
			{
				//set the texture coorinate
				glTexCoord2f(texturecoordinate[grid_cell_faces_in[i].texcoord[0] - 1]->u, texturecoordinate[grid_cell_faces_in[i].texcoord[0] - 1]->v);
			}
			glVertex3f(vertex[grid_cell_faces_in[i].faces[0] - 1]->x, vertex[grid_cell_faces_in[i].faces[0] - 1]->y, vertex[grid_cell_faces_in[i].faces[0] - 1]->z);


			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[grid_cell_faces_in[i].normal[1] - 1]->x, normals[grid_cell_faces_in[i].normal[1] - 1]->y, normals[grid_cell_faces_in[i].normal[1] - 1]->z);    //use them
			}
			if (istexture && materials.size() > 0 && materials[grid_cell_faces_in[i].mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[grid_cell_faces_in[i].texcoord[1] - 1]->u, texturecoordinate[grid_cell_faces_in[i].texcoord[1] - 1]->v);
			}
			glVertex3f(vertex[grid_cell_faces_in[i].faces[1] - 1]->x, vertex[grid_cell_faces_in[i].faces[1] - 1]->y, vertex[grid_cell_faces_in[i].faces[1] - 1]->z);


			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[grid_cell_faces_in[i].normal[2] - 1]->x, normals[grid_cell_faces_in[i].normal[2] - 1]->y, normals[grid_cell_faces_in[i].normal[2] - 1]->z);    //use them
			}
			if (istexture && materials.size() > 0 && materials[grid_cell_faces_in[i].mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[grid_cell_faces_in[i].texcoord[2] - 1]->u, texturecoordinate[grid_cell_faces_in[i].texcoord[2] - 1]->v);
			}
			glVertex3f(vertex[grid_cell_faces_in[i].faces[2] - 1]->x, vertex[grid_cell_faces_in[i].faces[2] - 1]->y, vertex[grid_cell_faces_in[i].faces[2] - 1]->z);


			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[grid_cell_faces_in[i].normal[3] - 1]->x, normals[grid_cell_faces_in[i].normal[3] - 1]->y, normals[grid_cell_faces_in[i].normal[3] - 1]->z);    //use them
			}
			if (istexture && materials.size() > 0 && materials[grid_cell_faces_in[i].mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[grid_cell_faces_in[i].texcoord[3] - 1]->u, texturecoordinate[grid_cell_faces_in[i].texcoord[3] - 1]->v);
			}
			glVertex3f(vertex[grid_cell_faces_in[i].faces[3] - 1]->x, vertex[grid_cell_faces_in[i].faces[3] - 1]->y, vertex[grid_cell_faces_in[i].faces[3] - 1]->z);
			
			glEnd();
		}
		else
		{
			glBegin(GL_TRIANGLES);
			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[grid_cell_faces_in[i].normal[0] - 1]->x, normals[grid_cell_faces_in[i].normal[0] - 1]->y, normals[grid_cell_faces_in[i].normal[0] - 1]->z);    //use them
			}
			if (istexture && materials.size() > 0 && materials[grid_cell_faces_in[i].mat]->texture != -1)  //if there are textures
			{
				//set the texture coorinate
				glTexCoord2f(texturecoordinate[grid_cell_faces_in[i].texcoord[0] - 1]->u, texturecoordinate[grid_cell_faces_in[i].texcoord[0] - 1]->v);
			}
			glVertex3f(vertex[grid_cell_faces_in[i].faces[0] - 1]->x, vertex[grid_cell_faces_in[i].faces[0] - 1]->y, vertex[grid_cell_faces_in[i].faces[0] - 1]->z);


			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[grid_cell_faces_in[i].normal[1] - 1]->x, normals[grid_cell_faces_in[i].normal[1] - 1]->y, normals[grid_cell_faces_in[i].normal[1] - 1]->z);    //use them
			}
			if (istexture && materials.size() > 0 && materials[grid_cell_faces_in[i].mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[grid_cell_faces_in[i].texcoord[1] - 1]->u, texturecoordinate[grid_cell_faces_in[i].texcoord[1] - 1]->v);
			}
			glVertex3f(vertex[grid_cell_faces_in[i].faces[1] - 1]->x, vertex[grid_cell_faces_in[i].faces[1] - 1]->y, vertex[grid_cell_faces_in[i].faces[1] - 1]->z);


			if (isnormals)   //if there are normals
			{
				glNormal3f(normals[grid_cell_faces_in[i].normal[2] - 1]->x, normals[grid_cell_faces_in[i].normal[2] - 1]->y, normals[grid_cell_faces_in[i].normal[2] - 1]->z);    //use them
			}
			if (istexture && materials.size() > 0 && materials[grid_cell_faces_in[i].mat]->texture != -1)
			{
				glTexCoord2f(texturecoordinate[grid_cell_faces_in[i].texcoord[2] - 1]->u, texturecoordinate[grid_cell_faces_in[i].texcoord[2] - 1]->v);
			}
			glVertex3f(vertex[grid_cell_faces_in[i].faces[2] - 1]->x, vertex[grid_cell_faces_in[i].faces[2] - 1]->y, vertex[grid_cell_faces_in[i].faces[2] - 1]->z);

		glEnd();
		}
	}
}

void objloader::DividFace2Grid()
{
	//对每个grid cell 中的面片的顶点做判断，有两个以上的顶点在cell中，就该面片属于cell
	for (int i = 0; i < divide_num; i++)
	{
		for (int j = 0; j < divide_num; j++)
		{
			GLfloat min_x = x_min + delt_width * j;
			GLfloat max_x = x_min + delt_width * (j+1);

			GLfloat min_z = z_min + delt_width * i;
			GLfloat max_z = z_min + delt_width * (i + 1);

			int vertex_in_grid = 0;
			vector<face> grid_face;
			for (int m = 0; m < faces.size(); m++)
			{
				if (faces[m]->four)
				{
					vertex_in_grid = 0;
					for (int n = 0; n < 4; n++)
					{
						if (vertex[faces[m]->faces[n]-1]->x >= min_x &&  vertex[faces[m]->faces[n]-1]->x <= max_x &&
							vertex[faces[m]->faces[n]-1]->z >= min_z &&  vertex[faces[m]->faces[n]-1]->z <= max_z)
						{
							vertex_in_grid++;
						}
					}
				}
				else
				{
					vertex_in_grid = 0;
					for (int n = 0; n < 3; n++)
					{
						if (vertex[faces[m]->faces[n] - 1]->x >= min_x &&  vertex[faces[m]->faces[n] - 1]->x <= max_x &&
							vertex[faces[m]->faces[n] - 1]->z >= min_z &&  vertex[faces[m]->faces[n] - 1]->z <= max_z)
						{
							vertex_in_grid++;
						}
					}
				}

				if (vertex_in_grid >= 1)
				{
					grid_face.push_back(*faces[m]);

					/*vector<face*>::iterator it = faces.begin() + m;
					faces.erase(it);*/
				}

			}//for (int m = 0; m < faces.size(); m++)

			if (grid_face.size() == 0)
			{
				grid_face.push_back(face(0, 0, 0, 0, 0, 0, 0, 0,0,0,0));
			}

			faces_in_grid.push_back(grid_face);

		}//for (int j = 0; j < divide_num; i++)
	}//for (int i = 0; i < divide_num; i++)

}

void objloader::clean()
{
	//delete all the dynamically allocated memory
	for (int i = 0; i < coord.size(); i++)
		delete coord[i];
	for (int i = 0; i < faces.size(); i++)
		delete faces[i];
	for (int i = 0; i < normals.size(); i++)
		delete normals[i];
	for (int i = 0; i < vertex.size(); i++)
		delete vertex[i];
	for (int i = 0; i < materials.size(); i++)
		delete materials[i];
	for (int i = 0; i < texturecoordinate.size(); i++)
		delete texturecoordinate[i];
	//and all elements from the vector
	coord.clear();
	faces.clear();
	normals.clear();
	vertex.clear();
	materials.clear();
	texturecoordinate.clear();
}

//load the filename textures (only BMP, R5G6B5 format)
unsigned int objloader::loadTexture(const char* filename)
{
	//nothing new in here
	unsigned int num;
	glGenTextures(1, &num);
	SDL_Surface* img = SDL_LoadBMP(filename);
	glBindTexture(GL_TEXTURE_2D, num);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img->w, img->h, 0, GL_RGB, GL_UNSIGNED_SHORT_5_6_5, img->pixels);
	glTexEnvi(GL_TEXTURE_2D, GL_TEXTURE_ENV_MODE, GL_MODULATE);       //maybe just this
	SDL_FreeSurface(img);
	texture.push_back(num);
	return num;
}

