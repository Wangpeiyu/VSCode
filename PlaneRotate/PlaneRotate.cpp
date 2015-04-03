// PlaneRotate.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <gl/glew.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <gl/freeglut.h>
#include <Eigen/Dense>
#include<stdio.h>
#include<math.h>
#include "ObjLoadView.h"

using namespace std;
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;

struct Vertex_3D
{
	GLfloat x;
	GLfloat y;
	GLfloat z;
};
/*向量：u=(u1,u2,u3) v=(v1,v2,v3)

叉积公式：u x v = { u2v3 - v2u3, u3v1 - v3u1, u1v2 - u2v1 }

点积公式：u * v = u1v1 + u2v2 + u3v33 = lul*lvl*COS(U, V)*/

GLfloat Vector_Dot_Multiple(const Vertex_3D &u, const Vertex_3D &v)
{
	return (u.x * v.x + u.y * v.y + u.z * v.z);
}

GLfloat Vector_Morl(const Vertex_3D &u)
{
	return (sqrt(u.x * u.x + u.y * u.y + u.z * u.z));
}

Vertex_3D Vector_Cross_Multiple(const Vertex_3D &u, const Vertex_3D &v)
{
	Vertex_3D result;
	result.x = u.y * v.z - u.z * v.y;
	result.y = u.z * v.x - u.x * v.z;
	result.z = u.x * v.y - u.y * v.x;

	return result;
}

void Vector3Normalize(Vertex_3D &vector_in)
{
	GLfloat length = sqrt(vector_in.x * vector_in.x + vector_in.y * vector_in.y + vector_in.z * vector_in.z);
	
	vector_in.x= vector_in.x / length;
	vector_in.y = vector_in.y / length;
	vector_in.z = vector_in.z / length;
}

void RotateArbitraryAxis(MatrixXf &M_Out, Vertex_3D &axis, const float &theta)
{
	Vector3Normalize(axis);
	float u = axis.x;
	float v = axis.y;
	float w = axis.z;

	M_Out(0,0) = cosf(theta) + (u * u) * (1 - cosf(theta));
	M_Out(0,1) = u * v * (1 - cosf(theta)) + w * sinf(theta);
	M_Out(0,2) = u * w * (1 - cosf(theta)) - v * sinf(theta);
	M_Out(0,3) = 0;

	M_Out(1,0) = u * v * (1 - cosf(theta)) - w * sinf(theta);
	M_Out(1,1) = cosf(theta) + v * v * (1 - cosf(theta));
	M_Out(1,2) = w * v * (1 - cosf(theta)) + u * sinf(theta);
	M_Out(1,3) = 0;

	M_Out(2,0) = u * w * (1 - cosf(theta)) + v * sinf(theta);
	M_Out(2,1) = v * w * (1 - cosf(theta)) - u * sinf(theta);
	M_Out(2,2) = cosf(theta) + w * w * (1 - cosf(theta));
	M_Out(2,3) = 0;

	M_Out(3,0) = 0;
	M_Out(3,1) = 0;
	M_Out(3,3) = 0;
	M_Out(3,3) = 1;
}

void RotateArbitraryLine(MatrixXf &M_Out, Vertex_3D v1, Vertex_3D v2, float theta)
{
	float a = v1.x;
	float b = v1.y;
	float c = v1.z;

	Vertex_3D p;
	p.x = v2.x - v1.x;
	p.y = v2.y - v1.y;
	p.z = v2.z - v1.z;

	Vector3Normalize(p);

	float u = p.x;
	float v = p.y;
	float w = p.z;

	float uu = u * u;
	float uv = u * v;
	float uw = u * w;
	float vv = v * v;
	float vw = v * w;
	float ww = w * w;
	float au = a * u;
	float av = a * v;
	float aw = a * w;
	float bu = b * u;
	float bv = b * v;
	float bw = b * w;
	float cu = c * u;
	float cv = c * v;
	float cw = c * w;

	float costheta = cosf(theta);
	float sintheta = sinf(theta);

	M_Out(0,0) = uu + (vv + ww) * costheta;
	M_Out(0,1) = uv * (1 - costheta) + w * sintheta;
	M_Out(0,2) = uw * (1 - costheta) - v * sintheta;
	M_Out(0,3) = 0;

	M_Out(1,0) = uv * (1 - costheta) - w * sintheta;
	M_Out(1,1) = vv + (uu + ww) * costheta;
	M_Out(1,2) = vw * (1 - costheta) + u * sintheta;
	M_Out(1,3) = 0;

	M_Out(2,0) = uw * (1 - costheta) + v * sintheta;
	M_Out(2,1) = vw * (1 - costheta) - u * sintheta;
	M_Out(2,2) = ww + (uu + vv) * costheta;
	M_Out(2,3) = 0;

	M_Out(3,0) = (a * (vv + ww) - u * (bv + cw)) * (1 - costheta) + (bw - cv) * sintheta;
	M_Out(3,1) = (b * (uu + ww) - v * (au + cw)) * (1 - costheta) + (cu - aw) * sintheta;
	M_Out(3,2) = (c * (uu + vv) - w * (au + bv)) * (1 - costheta) + (av - bu) * sintheta;
	M_Out(3,3) = 1;
}

Vertex_3D GetFaceNormal(Vertex_3D V1, Vertex_3D V2, Vertex_3D V3)
{
	Vertex_3D Normal;
	Normal.x = (V2.y - V1.y)*(V3.z - V1.z) - (V2.z - V1.z)*(V3.y - V1.y);
	Normal.y = (V3.x - V1.x)*(V2.z - V1.z) - (V2.x - V1.x)*(V3.z - V1.z);
	Normal.z = (V2.x - V1.x)*(V3.y - V1.y) - (V3.x - V1.x)*(V2.y - V1.y);

	return Normal;
}



struct point_t 
{
	double x, y;
};

struct circle_t 
{
	struct point_t center;
	double r;
};

int double_equals(double const a, double const b)
{
	static const double ZERO = 1e-9;
	return fabs(a - b) < ZERO;
}

double distance_sqr(struct point_t const* a, struct point_t const* b)
{
	return (a->x - b->x) * (a->x - b->x) + (a->y - b->y) * (a->y - b->y);
}

double distance(struct point_t const* a, struct point_t const* b)
{
	return sqrt(distance_sqr(a, b));
}

int insect(struct circle_t circles[], struct point_t points[])
{
	double d, a, b, c, p, q, r;
	double cos_value[2], sin_value[2];
	if (double_equals(circles[0].center.x, circles[1].center.x)
		&& double_equals(circles[0].center.y, circles[1].center.y)
		&& double_equals(circles[0].r, circles[1].r)) {
		return -1;
	}

	d = distance(&circles[0].center, &circles[1].center);
	if (d > circles[0].r + circles[1].r
		|| d < fabs(circles[0].r - circles[1].r)) {
		return 0;
	}

	a = 2.0 * circles[0].r * (circles[0].center.x - circles[1].center.x);
	b = 2.0 * circles[0].r * (circles[0].center.y - circles[1].center.y);
	c = circles[1].r * circles[1].r - circles[0].r * circles[0].r
		- distance_sqr(&circles[0].center, &circles[1].center);
	p = a * a + b * b;
	q = -2.0 * a * c;
	if (double_equals(d, circles[0].r + circles[1].r)
		|| double_equals(d, fabs(circles[0].r - circles[1].r))) {
		cos_value[0] = -q / p / 2.0;
		sin_value[0] = sqrt(1 - cos_value[0] * cos_value[0]);

		points[0].x = circles[0].r * cos_value[0] + circles[0].center.x;
		points[0].y = circles[0].r * sin_value[0] + circles[0].center.y;

		if (!double_equals(distance_sqr(&points[0], &circles[1].center),
			circles[1].r * circles[1].r)) {
			points[0].y = circles[0].center.y - circles[0].r * sin_value[0];
		}
		return 1;
	}

	r = c * c - b * b;
	cos_value[0] = (sqrt(q * q - 4.0 * p * r) - q) / p / 2.0;
	cos_value[1] = (-sqrt(q * q - 4.0 * p * r) - q) / p / 2.0;
	sin_value[0] = sqrt(1 - cos_value[0] * cos_value[0]);
	sin_value[1] = sqrt(1 - cos_value[1] * cos_value[1]);

	points[0].x = circles[0].r * cos_value[0] + circles[0].center.x;
	points[1].x = circles[0].r * cos_value[1] + circles[0].center.x;
	points[0].y = circles[0].r * sin_value[0] + circles[0].center.y;
	points[1].y = circles[0].r * sin_value[1] + circles[0].center.y;

	if (!double_equals(distance_sqr(&points[0], &circles[1].center),
		circles[1].r * circles[1].r)) {
		points[0].y = circles[0].center.y - circles[0].r * sin_value[0];
	}
	if (!double_equals(distance_sqr(&points[1], &circles[1].center),
		circles[1].r * circles[1].r)) {
		points[1].y = circles[0].center.y - circles[0].r * sin_value[1];
	}
	if (double_equals(points[0].y, points[1].y)
		&& double_equals(points[0].x, points[1].x)) {
		if (points[0].y > 0) {
			points[1].y = -points[1].y;
		}
		else {
			points[0].y = -points[0].y;
		}
	}
	return 2;
}

int main()
{
	struct circle_t circles[2];
	struct point_t points[2];

	std::cout << "请输入两圆x，y，半径(以逗号分开)：" << endl;

	std::cin >> circles[0].center.x >> circles[0].center.y >> circles[0].r >> circles[1].center.x >> circles[1].center.y >> circles[1].r;
	switch (insect(circles, points))
	{
	case -1:
		printf("THE CIRCLES ARE THE SAME/n");
		break;
	case 0:
		printf("NO INTERSECTION/n");
		break;
	case 1:
		printf("(%.3lf %.3lf)\n", points[0].x, points[0].y);
		break;
	case 2:
		printf("(%.3lf %.3lf) (%.3lf %.3lf)\n",
			points[0].x, points[0].y,
			points[1].x, points[1].y);
	}

	return 0;
}

int main(int argc, _TCHAR* argv[])
{
	objloader obj;
	obj.load("bunny_part.obj");

	struct circle_t circles[2];
	struct point_t points[2];

	Vertex_3D face_normal;
	Vertex_3D Plane_Normal;
	Vertex_3D Rotate_axis;
	Vertex_3D Rotate_axis_Normal;

	MatrixXf vertx_rotate_out(1,4);
	MatrixXf vertx_rotate_in(1,4);
	MatrixXf rotate_matrix(4,4);

	Plane_Normal.x = 0;
	Plane_Normal.y = 1;
	Plane_Normal.z = 0;

	for (int i = 0; i < obj.faces.size(); i++)
	{
		face_normal = GetFaceNormal(obj.vertex[obj.faces[]], vertexts[i * 3 + 1], vertexts[i * 3 + 2]);
	}

	
	Vector3Normalize(face_normal);

	Rotate_axis_Normal = Vector_Cross_Multiple(face_normal, Plane_Normal);
	float degree = acos(Vector_Dot_Multiple(face_normal, Plane_Normal) / (Vector_Morl(face_normal)* Vector_Morl(Plane_Normal)));

	Rotate_axis.x = vertexts[i * 3].x + Rotate_axis_Normal.x;
	Rotate_axis.y = vertexts[i * 3].y + Rotate_axis_Normal.y;
	Rotate_axis.z = vertexts[i * 3].z + Rotate_axis_Normal.z;


	RotateArbitraryLine(rotate_matrix, vertexts[i * 3], Rotate_axis, degree);

	cout << vertexts[i * 3].x << " " << vertexts[i * 3].y << "" << vertexts[i * 3].z << endl;
	vertx_rotate_in(0) = vertexts[i * 3 + j + 1].x;
	vertx_rotate_in(1) = vertexts[i * 3 + j + 1].y;
	vertx_rotate_in(2) = vertexts[i * 3 + j + 1].z;
	vertx_rotate_in(3) = 1;

	vertx_rotate_out = vertx_rotate_in * rotate_matrix;
	cout << vertx_rotate_out << endl;
	cout << endl;

	system("pause");
	return 0;
}

