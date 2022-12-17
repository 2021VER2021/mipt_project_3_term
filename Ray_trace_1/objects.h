#include "framework.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
//#include<amp.h>
//#include <omp.h>

// clock_1 used to measure time, while debuging
//#define clock_1
// OMP is the great CPU multithreading library, but 
// we found no improvement in computation speed 
//#define OMP
// GRAM you can define to do some crazy metric, but computation speed will be ruined
//#define GRAM

#ifdef OMP
int const MAX_THREADS = 3;
#endif

//  number of pixels, for each (axe x and y)
int const WIDTH = 400;
int const HEIGHT = 400;

// Window size parameters
int n1 = WIDTH;         
int n2 = HEIGHT;  

// max number of reflections, for each ray
int depth = 1;          

// How much pixels will be in one game pixel, per one side
int pixel = n1 / WIDTH; 

//params of FOVand scaling screen 
double  d = 1;                 
double  w = d;                
double h = (double)n2/n1 * d; 

// distance, after which we consider, that ray doesn't intersect any object
double const positive_inf = 10000000;

// parameter for intersection functions, need to define, when we intersect object
double epsilon = 0.00001;

// the distance, from which we will cast rays when calculating shadows
double ray_epsilon = 0.01;

// implementation of 3D vector
struct vector3 { 
    double x;
    double y;
    double z;
};

// symbol for real 3D space vector
using VEC3 = vector3; 

// math vector (eg in matrices)
using VEC = std::array<double, 3>; 

// data type, using to present each channel colour data
using COL_t = BYTE;

// container for RGB colour data
using COL = std::array<COL_t, 3>;

// class MATRIX using to rotate camera
using MATR = std::array<VEC, 3>;

#ifdef clock_1
clock_t s_1 = 0;
clock_t e_1 = 0;
#endif

// Gram matrix, we can use it, to change metric
MATR Gram = 
{
    std::array<double, 3>{1, 0, 0},
    std::array<double, 3>{0, 1, 0},
    std::array<double, 3>{0, 0, 1}
};

// colors in RGB format using container COL
const COL SADBROWN = {139, 69, 13};
const COL GREENYELLOW = {86, 127, 23};
const COL RED = { 250, 0, 0 };    
const COL BLUE = { 0, 0, 250 };
const COL YELLOW = { 250, 250, 0 };
const COL LIGHTCORAL = { 240, 128, 128 };
const COL TERRACOTBROWN = { 78, 59, 49 };

// if ray reachs positive_inf, we will color pixel in this color
COL BG_C = { 0, 191, 255 };

VEC3 operator - (VEC3 v1, VEC3 v2) {
    VEC3 V = { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    return V;
}

VEC3 operator - (VEC3 v1) {
    VEC3 V = { -v1.x, -v1.y, -v1.z };
    return V;
}

VEC3 operator + (VEC3 v1, VEC3 v2) {
    return v1 - (-v2);
}

VEC3 operator * (MATR m, VEC3 v1) {
    VEC3 v2 = { m[0][0] * v1.x + m[0][1] * v1.y + m[0][2] * v1.z,
               m[1][0] * v1.x + m[1][1] * v1.y + m[1][2] * v1.z,
               m[2][0] * v1.x + m[2][1] * v1.y + m[2][2] * v1.z};
    return v2;
}
MATR operator * (MATR m1, MATR m2) {
    MATR m3 = {
        std::array<double, 3>{
         m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0],
         m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1] + m1[0][2] * m2[2][1],
         m1[0][0] * m2[0][2] + m1[0][1] * m2[1][2] + m1[0][2] * m2[2][2]},
        {
            m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0] + m1[1][2] * m2[2][0],
            m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1] + m1[1][2] * m2[2][1],
            m1[1][0] * m2[0][2] + m1[1][1] * m2[1][2] + m1[1][2] * m2[2][2]
        },
        {
        m1[2][0] * m2[0][0] + m1[2][1] * m2[1][0] + m1[2][2] * m2[2][0],
        m1[2][0] * m2[0][1] + m1[2][1] * m2[1][1] + m1[2][2] * m2[2][1],
        m1[2][0] * m2[0][2] + m1[2][1] * m2[1][2] + m1[2][2] * m2[2][2]
        } 
    };
    return m3;
}

//#ifdef OMP
//double operator * (VEC3 v1, VEC3 v2) {
//    double answer = 0;
//#pragma omp parallel for reduction(+:answer)
//        for (int i = 0; i < 3; i++) {
//            answer += v1[i] * v2[i];
//        }
//#pragma omp barrier
//    return answer;
//}
//#else

double operator * (VEC3 v1, VEC3 v2) {
#ifdef GRAM
    v2 = Gram * v2;  // Gram matrix usage
#endif
    double answer = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    return answer;
}

//#endif

VEC3 operator * (double v1, VEC3 v2) {
    VEC3 V = { v2.x * v1, v2.y * v1, v2.z * v1 };
    return V;
}
VEC3 operator * (VEC3 v1, double v2) {
    return v2 * v1;
}
COL operator * (double c1, COL c2) {
    COL c = { static_cast<COL_t>(c2[0] * c1), static_cast<COL_t>(c2[1] * c1), static_cast<COL_t>(c2[2] * c1) };
    return c;
}
COL operator + (COL c1, COL c2) {
    return { static_cast<COL_t>(c1[0] + c2[0]), static_cast<COL_t>(c1[1] + c2[1]), static_cast<COL_t>(c1[2] + c2[2]) };
}

/// <summary>
/// 
/// </summary>
/// <param name="v1 - input vector"></param>
/// <returns> length of vector </returns>
double abs(VEC3 &v1) {
    double V = v1 * v1;
    return std::sqrt(V);
}

/// <summary>
/// 
/// </summary>
/// <param name="v1 - first vector"></param>
/// <param name="v2 - second vector"></param>
/// <returns> cross product of v1 and v2 </returns>
VEC3 cross(VEC3 v1, VEC3 v2) {
    VEC3 V = { v1.y * v2.z - v1.z * v2.y, v2.x * v1.z - v1.x * v2.z, v1.x * v2.y - v2.x * v1.y };
    return V;
}

/// <summary>
/// 
/// </summary>
/// <param name="v"> - input vector, length != 0 </param>
/// <returns> normalized vector, length = 1</returns>
VEC3 normalize(VEC3& v) {
    return (1 / abs(v)) * v;
}

/// <summary>
/// 
/// </summary>
/// <param name="v"> - input vector, length != 0 </param>
/// <returns> normalized vector, length = 1</returns>
VEC3 normalize(VEC3&& v) {
    VEC3 v1 = (1 / abs(v)) * v;
    return v1;
}

// Important part of orientation, vector pointing up
VEC3 UP = normalize({ 0, 1, 0 });

/// <summary>
///  Rotation around x axe
/// </summary>
/// <param name="angle, radians"></param>
/// <param name="v - input vector"></param>
/// <returns> vector v, but rotated by angle radian around x </returns>
VEC3 xRotate(double angle, VEC3 v) {   // FIXIT
    MATR xRotation = { 
        VEC{1, 0, 0},
        VEC{0, std::cos(angle), -std::sin(angle)},
        VEC{ 0, std::sin(angle), std::cos(angle)}};
    return xRotation * v;
}

/// <summary>
///  Rotation around y axe
/// </summary>
/// <param name="angle, radians"></param>
/// <param name="v - input vector"></param>
/// <returns> vector v, but rotated by angle radianm around y </returns>
VEC3 yRotate(double angle, VEC3 v) {    // FIXIT
    MATR yRotation = { 
        VEC{std::cos(angle), 0, std::sin(angle)},
        VEC{0, 1, 0},
        VEC{ -std::sin(angle), 0, std::cos(angle)}};
    return yRotation * v;
}

/// <summary>
///  Rotation around z axe
/// </summary>
/// <param name="angle, radians"></param>
/// <param name="v - input vector"></param>
/// <returns> vector v, but rotated by angle radian around z </returns>
VEC3 zRotate(double angle, VEC3 v) {    // FIXIT
    MATR zRotation = {
        VEC{std::cos(angle), -std::sin(angle), 0},
        VEC{std::sin(angle), std::cos(angle), 0},
        VEC{0, 0, 1} };
    return zRotation * v;
}

/// <summary>
///  Rotation around x axe, and y axe
/// </summary>
/// <param name="angle, radians"></param>
/// <param name="v - input vector"></param>
/// <returns> vector v, but rotated by angle_x and angle_y around x and y </returns>
VEC3 xyRotate(double angle_x, double angle_y, VEC3 v) { // FIXIT
    MATR xRotation = {
        VEC{1, 0, 0},
        VEC{0, std::cos(angle_x), -std::sin(angle_x)},
        VEC{ 0, std::sin(angle_x), std::cos(angle_x)} };
    MATR yRotation = { 
        VEC{std::cos(angle_y), 0, std::sin(angle_y)},
        VEC{0, 1, 0},
        VEC{ -std::sin(angle_y), 0, std::cos(angle_y)} };
    return xRotation * yRotation * v;
}

/// <summary>
///  Helpful function, to calculate reflections
/// </summary>
/// <param name="R - falling ray"></param>
/// <param name="N - normal vector to surface"></param>
/// <returns> reflected ray </returns>
VEC3 ReflectRay(VEC3& R, VEC3& N) {
    return  2 * (R * N) * N - R;
}

/// <summary>
///  Helpful function, to calculate reflections
/// </summary>
/// <param name="R - falling ray"></param>
/// <param name="N - normal vector to surface"></param>
/// <returns> reflected ray </returns>
VEC3 ReflectRay(VEC3&& R, VEC3& N) {
    return  2 * (R * N) * N - R;
}

/// <summary>
/// Actualy, this is a kind of magic, why
/// this strange method of setting rectangles on 
/// canvas is the best in performance
/// </summary>
/// <param name="hdc - something, that allow you to use graphics"></param>
/// <param name="rect - reference to the RECT, wich you want to paint"></param>
/// <param name="colour - RGB color of the rect"></param>
void PaintRect(HDC hdc, RECT* rect, COLORREF colour)
{
    COLORREF oldcr = SetBkColor(hdc, colour);
    ExtTextOut(hdc, 0, 0, ETO_OPAQUE, rect, (LPCWSTR)"", 0, 0);
    SetBkColor(hdc, oldcr);
}

/// FIXIT?
/// <summary>
/// This class for discribe light
/// RAII done, but Inheritance...
/// </summary>
class LightObj {
protected:
    VEC3* d;           // direction (for direct) or position (for point)
    char* type;
    double* intensity;
public:
    /// <summary>
    /// basic constructors (brand new data)
    /// </summary>
    /// <param name="v"></param>
    /// <param name="s"></param>
    /// <param name="i"></param>
    LightObj(VEC3 const& v, char const& s, double const& i) {
        d = new VEC3; *d = v;
        type = new char; *type = s;
        intensity = new double; *intensity = i;
    };
    LightObj(LightObj const& obj) : LightObj(*obj.d, *obj.type, *obj.intensity) {};
    LightObj(char s, double i) : LightObj({ 0, 0, 0 }, s, i) {};

    /// <summary>  
    ///  Right value constructor
    /// </summary>
    /// <param name="t"></param>
    LightObj(LightObj&& t) {
        d = t.d;                 t.d = nullptr;
        type = t.type;           t.type = nullptr;
        intensity = t.intensity; t.intensity = nullptr;
    }

    ~LightObj() {
        delete d;
        delete type;
        delete intensity;
    }

    /// <summary>
    /// Copy assignment
    /// </summary>
    /// <param name="s"></param>
    /// <returns></returns>
    LightObj& operator = (LightObj const& s) {
        LightObj tmp(s);
        std::swap(this->d, tmp.d);
        std::swap(this->type, tmp.type);
        std::swap(this->intensity, tmp.intensity);
        return *this;
    }
    /// <summary>
    /// moving assignment
    /// </summary>
    /// <param name="t"></param>
    /// <returns></returns>
    LightObj& operator = (LightObj&& t) {
        LightObj tmp(std::move(t));
        std::swap(this->d, t.d);
        std::swap(this->type, t.type);
        std::swap(this->intensity, t.intensity);
        return *this;
    }
    VEC3 get_d() {
        return *d;
    }
    char get_type() {
        return *type;
    }
    double get_intensity() {
        return *intensity;
    }
};
/// <summary>
/// GenericObject Interface for every object you can visualise 
/// RAII done
/// </summary>
class GenericObject {
protected:
    COL color;

    // how it matt or glossy (-1 - matt, 0-10000 - glossy)
    double reflective; 

    // reflections on/off and how effective
    double specular;  
public:
    GenericObject(COL& c, double& refl, double& spec) {
        color = c;
        reflective = refl;
        specular = spec;
    };
    GenericObject(GenericObject& t) : GenericObject(t.color, t.reflective, t.specular) {}
    GenericObject(GenericObject&& t) {
        color = t.color;
        reflective = t.reflective;
        specular = t.specular;
    }
    ~GenericObject() {
    }
    /// <summary>
    /// Function calculates closest intercection for ray (actually line) and object 
    /// </summary>
    /// <param name="O - camera center in affine space"></param> 
    /// <param name="V - vector for the Intercecting ray"></param> 
    /// <param name="t_min - minimum distance for intercection"></param>
    /// <param name="t_max - maximum distance for intercection"></param>
    /// <param name="intersection - function output variable"></param>
    virtual void Intercections(VEC3& O, VEC3& V, double& t_min, double& t_max, double& intersection) = 0;
    virtual VEC3 get_norm(VEC3 P) = 0;
    virtual double get_reflective() = 0;
    virtual double get_specular() = 0;
    virtual COL get_color() = 0;
};

/// <summary>
/// Object, that discribe spheres 
/// RAII done
/// </summary>
class SphereObj : public GenericObject {
protected:
    VEC3 center;
    double radius;
public:
    SphereObj(VEC3 v, double r, COL c, double spec, double refl) : GenericObject(c, refl, spec) {
        radius = r;
        center = v;
    };
    SphereObj(VEC3 v, double r, COL c, double spec) : SphereObj(v, r, c, spec, 0.0) {};
    SphereObj(VEC3 v, double r, COL c) : SphereObj(v, r, c, -1.0) {};

    SphereObj(SphereObj const& t) : SphereObj(t.center, t.radius, t.color, t.specular, t.reflective) {};

    SphereObj(SphereObj&& t) : GenericObject(t) {// may not work as intended
        center = t.center;
        radius = t.radius;
    }

    ~SphereObj() {
    };


    /// <summary>
    /// Copy assignment
    /// </summary>
    /// <param name="s"></param>
    /// <returns></returns>
    SphereObj& operator = (SphereObj const& s) {
        SphereObj tmp(s);
        std::swap(this->center, tmp.center);
        std::swap(this->radius, tmp.radius);
        std::swap(this->color, tmp.color);
        std::swap(this->reflective, tmp.reflective);
        std::swap(this->specular, tmp.specular);
        return *this;
    }

    /// <summary>
    /// moving assignment
    /// </summary>
    /// <param name="t"></param>
    /// <returns></returns>
    SphereObj& operator = (SphereObj&& t) {
        SphereObj tmp(std::move(t));
        std::swap(this->center, t.center);
        std::swap(this->radius, t.radius);
        std::swap(this->color, t.color);
        std::swap(this->reflective, t.reflective);
        std::swap(this->specular, t.specular);
        return *this;
    }

    virtual void Intercections(VEC3& O, VEC3& V, double& t_min, double& t_max, double& closest_t) override {
        VEC3 C = center;
        double r = radius;
        VEC3 OC = O - C;
        double k1 = V * V;
        double k2 = 2 * (OC * V);
        double k3 = (OC * OC) - r * r;
        double diskr = k2 * k2 - 4 * k1 * k3;

        closest_t = positive_inf;
        double t1, t2;

        if (diskr < 0) {
            t1 = positive_inf;
            t2 = positive_inf;
        }
        else {
            t1 = (-k2 + std::sqrt(diskr)) / (2 * k1);
            t2 = (-k2 - std::sqrt(diskr)) / (2 * k1);
        }
        if ((t1 >= t_min) && (t1 <= t_max)) {
            closest_t = t1;
        }
        if ((t2 >= t_min) && (t2 <= t_max)) {
            closest_t = t2;
        }
    };

    virtual COL get_color() override {
        return GenericObject::color;
    }
    virtual double get_reflective() override {
        return GenericObject::reflective;
    }
    virtual double get_specular() override {
        return GenericObject::specular;
    }
    VEC3 get_center() {
        return center;
    }
    void set_center(VEC3 p) {
        center = p;
    }
    virtual VEC3 get_norm(VEC3 P) override {
        VEC3 N = P - center;
        return N;
    }
};

/// <summary>
/// Object, that discribe planes
/// RAII done
/// </summary>
class PlaneObj : public GenericObject {
protected:
    VEC3 norm;

    // the point, that belong the plane
    VEC3 param;
public:
    PlaneObj(VEC3 norm, VEC3 param, COL c, double spec, double refl) : GenericObject(c, refl, spec) {
        PlaneObj::norm = norm;
        PlaneObj::param = param;
    }
    PlaneObj(PlaneObj const& t) : PlaneObj(t.norm, t.param, t.color, t.specular, t.reflective) {};
    PlaneObj(PlaneObj&& t) : GenericObject(t) {
        norm = t.norm;
        param = t.param;
    }
    ~PlaneObj() {
    }

    /// <summary>
    /// Copy assignment
    /// </summary>
    /// <param name="s"></param>
    /// <returns></returns>
    PlaneObj& operator = (PlaneObj const& s) {
        PlaneObj tmp(s);
        std::swap(this->norm, tmp.norm);
        std::swap(this->param, tmp.param);
        std::swap(this->color, tmp.color);
        std::swap(this->reflective, tmp.reflective);
        std::swap(this->specular, tmp.specular);
        return *this;
    }
    /// <summary>
    /// moving assignment
    /// </summary>
    /// <param name="t"></param>
    /// <returns></returns>
    PlaneObj& operator = (PlaneObj&& t) {
        PlaneObj tmp(std::move(t));
        std::swap(this->param, t.param);
        std::swap(this->norm, t.norm);
        std::swap(this->color, t.color);
        std::swap(this->reflective, t.reflective);
        std::swap(this->specular, t.specular);
        return *this;
    }

    virtual void Intercections(VEC3& O, VEC3& V, double& t_min, double& t_max, double& closest_t) override {
        double dot = norm * V;
        double t1 = positive_inf;
        closest_t = positive_inf;
        if (abs(dot) >= epsilon) {
            VEC3 W = O-param;
            double fac = -(norm * W)/dot;
            if (fac > 0) {
                t1 = fac;
            }
        }
        if ((t1 >= t_min) && (t1 <= t_max)) {
            closest_t = t1;
        }
    }

    virtual VEC3 get_norm(VEC3 P) override {
        return norm;
    }
    virtual COL get_color() override {
        return GenericObject::color;
    }
    virtual double get_reflective() override {
        return GenericObject::reflective;
    }
    virtual double get_specular() override {
        return GenericObject::specular;
    }
};

/// <summary>
/// RAII DONE
/// </summary>
class TriangleObj : public GenericObject {
protected:
    VEC3 norm;
    VEC3 point_1;  //FIXED
    VEC3 point_2;
    VEC3 point_3;
public:
    TriangleObj(VEC3 Point_1, VEC3 Point_2, VEC3 Point_3, COL c, double spec, double refl): GenericObject(c, refl, spec) {
        norm = cross(Point_2 - Point_1, Point_3 - Point_1);
        point_1 = Point_1;
        point_2 = Point_2;
        point_3 = Point_3;
    }
    TriangleObj(TriangleObj const& t) : TriangleObj(t.point_1, t.point_2, t.point_3, t.color, t.specular, t.reflective) {};
    TriangleObj(TriangleObj&& t) : GenericObject(t) {
        point_1 = t.point_1;
        point_2 = t.point_2;
        point_2 = t.point_3;
        norm = t.norm;
    }
    ~TriangleObj() {
    }

    /// <summary>
    /// Copy assignment
    /// </summary>
    /// <param name="s"></param>
    /// <returns></returns>
    TriangleObj& operator = (TriangleObj const& s) {
        TriangleObj tmp(s);
        std::swap(this->norm, tmp.norm);
        std::swap(this->point_1, tmp.point_1);
        std::swap(this->point_2, tmp.point_2);
        std::swap(this->point_3, tmp.point_3);
        std::swap(this->color, tmp.color);
        std::swap(this->reflective, tmp.reflective);
        std::swap(this->specular, tmp.specular);
        return *this;
    }

    /// <summary>
    /// moving assignment
    /// </summary>
    /// <param name="t"></param>
    /// <returns></returns>
    TriangleObj& operator = (TriangleObj&& t) {
        TriangleObj tmp(std::move(t));
        std::swap(this->point_1, t.point_1);
        std::swap(this->point_2, t.point_2);
        std::swap(this->point_3, t.point_3);
        std::swap(this->norm, t.norm);
        std::swap(this->color, t.color);
        std::swap(this->reflective, t.reflective);
        std::swap(this->specular, t.specular);
        return *this;
    }

    virtual void Intercections(VEC3& O, VEC3& V, double& t_min, double& t_max, double& closest_t) override {
        VEC3 e1 = point_2 - point_1;
        VEC3 e2 = point_3 - point_1;
        VEC3 pvec = cross(V, e2);
        double det = (e1 * pvec);
        double t = positive_inf;
        closest_t = positive_inf;
        if (!(det < epsilon && det > -epsilon)) {
            double inv_det = 1 / det;
            VEC3 tvec = O - point_1;
            double u = (tvec * pvec) * inv_det;
            if (!(u < 0 || u > 1)) {
                VEC3 qvec = cross(tvec, e1);
                double v = (V * qvec) * inv_det;
                if (!(v < 0 || v + u > 1)) {
                    t = (e2 * qvec) * inv_det;
                }
            }
        }
        if ((t >= t_min) && (t <= t_max)) {
            closest_t = t;
        }
    }

    virtual VEC3 get_norm(VEC3 P) override {
        return norm;
    }
    virtual COL get_color() override {
        return GenericObject::color;
    }
    virtual double get_reflective() override {
        return GenericObject::reflective;
    }
    virtual double get_specular() override {
        return GenericObject::specular;
    }
};
/// <summary>
///  RAII DONE
/// </summary>
class RectangleObj : public GenericObject {
protected:
    VEC3 point_0;
    VEC3 point_1;  
    VEC3 point_2;
    VEC3 point_3;
    VEC3 norm;
public:
    /// <summary>
    /// point format:
    /// 1-----2
    /// |-----|
    /// 3-----4
    /// </summary>
    RectangleObj(VEC3 Point_1, VEC3 Point_2, VEC3 Point_3, VEC3 Point_4, COL c, double spec, double refl) : GenericObject(c, refl, spec) {
        point_0 = Point_1;
        point_1 = Point_2;
        point_2 = Point_3;
        point_3 = Point_4;
        norm = normalize(cross(Point_2 - Point_1, Point_3 - Point_1));
    }
    RectangleObj(RectangleObj const& t) : RectangleObj(t.point_0, t.point_1, t.point_2, t.point_3, t.color, t.specular, t.reflective) {};
    RectangleObj(RectangleObj&& t) : RectangleObj(t) {
        point_0 = t.point_0;
        point_1 = t.point_1;
        point_2 = t.point_2;
        point_3 = t.point_3;
        norm = t.norm;
    }
    ~RectangleObj() {}

    /// <summary>
    /// Copy assignment
    /// </summary>
    /// <param name="s"></param>
    /// <returns></returns>
    RectangleObj& operator = (RectangleObj const& s) {
        RectangleObj tmp(s);
        std::swap(this->norm, tmp.norm);
        std::swap(this->point_0, tmp.point_0);
        std::swap(this->point_1, tmp.point_1);
        std::swap(this->point_2, tmp.point_2);
        std::swap(this->point_3, tmp.point_3);
        std::swap(this->color, tmp.color);
        std::swap(this->reflective, tmp.reflective);
        std::swap(this->specular, tmp.specular);
        return *this;
    }

    /// <summary>
    /// moving assignment
    /// </summary>
    /// <param name="t"></param>
    /// <returns></returns>
    RectangleObj& operator = (RectangleObj&& t) {
        RectangleObj tmp(std::move(t));
        std::swap(this->point_0, t.point_0);
        std::swap(this->point_1, t.point_1);
        std::swap(this->point_2, t.point_2);
        std::swap(this->point_3, t.point_3);
        std::swap(this->norm, t.norm);
        std::swap(this->color, t.color);
        std::swap(this->reflective, t.reflective);
        std::swap(this->specular, t.specular);
        return *this;
    }

    void TriIntersections(VEC3& O, VEC3& V, double& t_min, double& t_max, VEC3&& e1, VEC3&& e2, VEC3& point, double& closest_t) {
        VEC3 pvec = cross(V, e2);
        double det = (e1 * pvec);
        double t = positive_inf;
        if (!(det < epsilon && det > -epsilon)) {
            double inv_det = 1 / det;
            VEC3 tvec = O - point;
            double u = (tvec * pvec) * inv_det;
            if (!(u < 0 || u > 1)) {
                VEC3 qvec = cross(tvec, e1);
                double v = (V * qvec) * inv_det;
                if (!(v < 0 || v + u > 1)) {
                    t = (e2 * qvec) * inv_det;
                }
            }
        }
        if (((t >= t_min) && (t <= t_max)) && (t < closest_t)) {
            closest_t = t;
        }
    }

    virtual void Intercections(VEC3& O, VEC3& V, double& t_min, double& t_max, double& closest_t) override {
        closest_t = positive_inf;
        TriIntersections(O, V, t_min, t_max, point_1 - point_0, point_2 - point_0, point_0, closest_t);
        TriIntersections(O, V, t_min, t_max, point_2 - point_1, point_3 - point_1, point_1, closest_t);
    }

    virtual VEC3 get_norm(VEC3 P) override {
        return -norm;
    }
    virtual COL get_color() override {
        return GenericObject::color;
    }
    virtual double get_reflective() override {
        return GenericObject::reflective;
    }
    virtual double get_specular() override {
        return GenericObject::specular;
    }
};
/// <summary>
/// point_format:
/// -----2
/// |    |
/// 1-----
/// RAII done
/// </summary>
class WallObj : public RectangleObj {
public:
    WallObj(VEC3 Point_1, VEC3 Point_2, COL c, double spec, double refl) : RectangleObj(
        Point_1,
        Point_1 + ((Point_2 - Point_1) - (UP * (Point_2 - Point_1)) * UP),
        Point_1 + (UP * (Point_2 - Point_1)) * UP,
        Point_2,
        c, spec, refl) {}
    WallObj(WallObj const& t) : RectangleObj(t.point_0, t.point_1, t.point_2, t.point_3, t.color, t.specular, t.reflective) {};
 
   WallObj(WallObj&& t) : WallObj(t) {
        point_0 = t.point_0;
        point_1 = t.point_1;
        point_2 = t.point_2;
        point_3 = t.point_3;
        norm = t.norm;
    }

    /// <summary>
    /// Copy assignment
    /// </summary>
    /// <param name="s"></param>
    /// <returns></returns>
    WallObj& operator = (WallObj const& s) {
        WallObj tmp(s);
        std::swap(this->norm, tmp.norm);
        std::swap(this->point_0, tmp.point_0);
        std::swap(this->point_1, tmp.point_1);
        std::swap(this->point_2, tmp.point_2);
        std::swap(this->point_3, tmp.point_3);
        std::swap(this->color, tmp.color);
        std::swap(this->reflective, tmp.reflective);
        std::swap(this->specular, tmp.specular);
        return *this;
    }

    /// <summary>
    /// moving assignment
    /// </summary>
    /// <param name="t"></param>
    /// <returns></returns>
    WallObj& operator = (WallObj&& t) {
        WallObj tmp(std::move(t));
        std::swap(this->point_0, t.point_0);
        std::swap(this->point_1, t.point_1);
        std::swap(this->point_2, t.point_2);
        std::swap(this->point_3, t.point_3);
        std::swap(this->norm, t.norm);
        std::swap(this->color, t.color);
        std::swap(this->reflective, t.reflective);
        std::swap(this->specular, t.specular);
        return *this;
    }
};
