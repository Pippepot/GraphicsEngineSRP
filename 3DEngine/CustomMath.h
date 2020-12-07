#pragma once

#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <math.h>
using namespace std;

struct mat4x4
{

	float m[4][4] = { 0 };

	mat4x4& operator * (mat4x4& m2)
	{
		mat4x4 matrix;
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				matrix.m[r][c] = m[r][0] * m2.m[0][c] + m[r][1] * m2.m[1][c] + m[r][2] * m2.m[2][c] + m[r][3] * m2.m[3][c];
		return matrix;
	}

};

struct vec3d
{
public:

	float x;
	float y;
	float z;
	float w; // Need a 4th term to perform sensible matrix vector multiplication

	vec3d()
	{
		x = y = z = 0;
		w = 1;
	}

	vec3d(float a, float b, float c, float d = 1)
	{
		x = a;
		y = b;
		z = c;
		w = d;
	}

	float Length() const
	{
		return sqrtf(x * x + y * y + z * z);
	}

	vec3d Normalise()
	{
		return *this / Length();
	}

#pragma region Operators


	vec3d& operator+=(const vec3d& rhs)
	{
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	};

	vec3d& operator-=(const vec3d& rhs)
	{
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		return *this;
	};

	vec3d& operator*=(const vec3d& rhs)
	{
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
		return *this;
	};

	vec3d& operator/=(const vec3d& rhs)
	{
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
		return *this;
	};

	vec3d& operator+(const vec3d& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x + rhs.x;
		r.y = this->y + rhs.y;
		r.z = this->z + rhs.z;
		return r;
	};

	vec3d& operator-(const vec3d& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x - rhs.x;
		r.y = this->y - rhs.y;
		r.z = this->z - rhs.z;
		return r;
	};

	vec3d& operator*(const vec3d& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x * rhs.x;
		r.y = this->y * rhs.y;
		r.z = this->z * rhs.z;
		return r;
	};

	vec3d& operator/(const vec3d& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x / rhs.x;
		r.y = this->y / rhs.y;
		r.z = this->z / rhs.z;
		return r;
	};


	vec3d& operator+=(const float& rhs)
	{
		this->x += rhs;
		this->y += rhs;
		this->z += rhs;
		return *this;
	};

	vec3d& operator-=(const float& rhs)
	{
		this->x -= rhs;
		this->y -= rhs;
		this->z -= rhs;
		return *this;
	};

	vec3d& operator*=(const float& rhs)
	{
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	};

	vec3d& operator/=(const float& rhs)
	{
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
		return *this;
	}

	vec3d& operator+(const float& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x + rhs;
		r.y = this->y + rhs;
		r.z = this->z + rhs;
		return r;
	};

	vec3d& operator-(const float& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x - rhs;
		r.y = this->y - rhs;
		r.z = this->z - rhs;
		return r;
	};

	vec3d& operator*(const float& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x * rhs;
		r.y = this->y * rhs;
		r.z = this->z * rhs;
		return r;
	};

	vec3d& operator/(const float& rhs)
	{
		vec3d r = { 0,0,0 };
		r.x = this->x / rhs;
		r.y = this->y / rhs;
		r.z = this->z / rhs;
		return r;
	};



#pragma endregion



};

struct vec2d
{
	float u = 0;
	float v = 0;
	float w = 1;

	vec2d()
	{
		u = v = 0;
		w = 1;
	}

	vec2d(float a, float b, float c = 1)
	{
		u = a;
		v = b;
		w = c;
	}

#pragma region Operators


	vec2d& operator+=(const vec2d& rhs)
	{
		this->u += rhs.u;
		this->v += rhs.v;
		return *this;
	};

	vec2d& operator-=(const vec2d& rhs)
	{
		this->u -= rhs.u;
		this->v -= rhs.v;
		return *this;
	};

	vec2d& operator*=(const vec2d& rhs)
	{
		this->u *= rhs.u;
		this->v *= rhs.v;
		return *this;
	};

	vec2d& operator/=(const vec2d& rhs)
	{
		this->u *= rhs.u;
		this->v *= rhs.v;
		return *this;
	};

	vec2d& operator+(const vec2d& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u + rhs.u;
		r.v = this->v + rhs.v;
		return r;
	};

	vec2d& operator-(const vec2d& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u - rhs.u;
		r.v = this->v - rhs.v;
		return r;
	};

	vec2d& operator*(const vec2d& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u * rhs.u;
		r.v = this->v * rhs.v;
		return r;
	};

	vec2d& operator/(const vec2d& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u / rhs.u;
		r.v = this->v / rhs.v;
		return r;
	};


	vec2d& operator+=(const float& rhs)
	{
		this->u += rhs;
		this->v += rhs;
		return *this;
	};

	vec2d& operator-=(const float& rhs)
	{
		this->u -= rhs;
		this->v -= rhs;
		return *this;
	};

	vec2d& operator*=(const float& rhs)
	{
		this->u *= rhs;
		this->v *= rhs;
		return *this;
	};

	vec2d& operator/=(const float& rhs)
	{
		this->u /= rhs;
		this->v /= rhs;
		return *this;
	}

	vec2d& operator+(const float& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u + rhs;
		r.v = this->v + rhs;
		return r;
	};

	vec2d& operator-(const float& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u - rhs;
		r.v = this->v - rhs;
		return r;
	};

	vec2d& operator*(const float& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u * rhs;
		r.v = this->v * rhs;
		return r;
	};

	vec2d& operator/(const float& rhs)
	{
		vec2d r = { 0,0 };
		r.u = this->u / rhs;
		r.v = this->v / rhs;
		return r;
	};



#pragma endregion



};

struct triangle
{
	vec3d p[3];
	vec2d t[3];

	olc::Pixel col;
};

float Vector_Dot(vec3d& v1, vec3d& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

vec3d Vector_Cross(vec3d& v1, vec3d& v2)
{
	vec3d v;
	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;
	return v;
}

vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd, float &t)
{
	plane_n = plane_n.Normalise();
	float plane_d = -Vector_Dot(plane_n, plane_p);
	float ad = Vector_Dot(lineStart, plane_n);
	float bd = Vector_Dot(lineEnd, plane_n);
	t = (-plane_d - ad) / (bd - ad);
	vec3d lineStartToEnd = lineEnd - lineStart;
	vec3d lineToIntersect = lineStartToEnd * t;
	return lineStart + lineToIntersect;
}



vec3d Matrix_MultiplyVector(mat4x4& m, vec3d& i)
{
	vec3d v;
	v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
	v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
	v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
	v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
	return v;
}


mat4x4 Matrix_Identity()
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
};

mat4x4 Matrix_RotationX(float fAngleRad)
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[1][2] = sinf(fAngleRad);
	matrix.m[2][1] = -sinf(fAngleRad);
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0f;
	return matrix;
};

mat4x4 Matrix_RotationY(float fAngleRad)
{
	mat4x4 matrix;
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][2] = sinf(fAngleRad);
	matrix.m[2][0] = -sinf(fAngleRad);
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0f;
	return matrix;
};

mat4x4 Matrix_RotationZ(float fAngleRad)
{
	mat4x4 matrix;
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][1] = sinf(fAngleRad);
	matrix.m[1][0] = -sinf(fAngleRad);
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
};

mat4x4 Matrix_Translation(float x, float y, float z)
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	matrix.m[3][0] = x;
	matrix.m[3][1] = y;
	matrix.m[3][2] = z;
	return matrix;
};

mat4x4 Matrix_Projection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
{
	float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
	mat4x4 matrix;
	matrix.m[0][0] = fAspectRatio * fFovRad;
	matrix.m[1][1] = fFovRad;
	matrix.m[2][2] = fFar / (fFar - fNear);
	matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matrix.m[2][3] = 1.0f;
	matrix.m[3][3] = 0.0f;
	return matrix;
};

mat4x4 Matrix_PointAt(vec3d& pos, vec3d& target, vec3d& up)
{
	// Calculate forward direction
	vec3d newForward = target - pos;
	newForward = newForward.Normalise();

	// Calculate new up direction
	vec3d a = newForward * Vector_Dot(up, newForward);
	vec3d newUp = up - a;
	newUp = newUp.Normalise();

	// Calculate new right direction
	vec3d newRight = Vector_Cross(newUp, newForward);

	// Construct Dimensioning and Translation Matrix	
	mat4x4 matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;

}

mat4x4 Matrix_QuickInverse(mat4x4& m) // Only for Rotation/Translation Matrices
{
	mat4x4 matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
	matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
	matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
	matrix.m[3][3] = 1.0f;
	return matrix;
}
