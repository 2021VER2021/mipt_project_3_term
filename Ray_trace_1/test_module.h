#include <iostream>
#include <vector>
#include <cmath>
#include <array>

int n1 = 300;  //  rename?
int n2 = 300;
int pixel = 8;
double  d = 1;  // params of FOV and scaling screen // how does it work actually??
double  w = d;  // params of FOV and scaling screen
double h = (double)n2/n1 * d;  // params of FOV and scaling screen
double const positive_inf = 100000000;
double epsilon = 0.00001;
using VEC = std::array<double, 3>;
using COL_t = BYTE;
using COL = std::array<COL_t, 3>;
using MATR = std::array<VEC, 3>;

/// <summary>
/// Gram moment
/// </summary>
MATR Gram = {std::array<double, 3>{1, 0, 0}, {0, 1, 0 }, { 0, 0, 1 }};

const COL RED = { 250, 0, 0 };    // Map? maybe
const COL BLUE = { 0, 0, 250 };
const COL YELLOW = { 250, 250, 0 };
const COL BG_C = { 0, 191, 255 };

VEC operator - (VEC v1, VEC v2) {
    VEC V = { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
    return V;
}
VEC operator - (VEC v1) {
    VEC V = { -v1[0], -v1[1], -v1[2] };
    return V;
}
VEC operator + (VEC v1, VEC v2) {
    return v1 - (-v2);
}

VEC operator * (MATR m, VEC v1) {
    VEC v2 = { m[0][0] * v1[0] + m[0][1] * v1[1] + m[0][2] * v1[2],
               m[1][0] * v1[0] + m[1][1] * v1[1] + m[1][2] * v1[2],
               m[2][0] * v1[0] + m[2][1] * v1[1] + m[2][2] * v1[2]};
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

double operator * (VEC v1, VEC v2) {
    v2 = Gram * v2;
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
VEC operator * (double v1, VEC v2) {
    VEC V = { v2[0] * v1, v2[1] * v1, v2[2] * v1 };
    return V;
}
VEC operator * (VEC v1, double v2) {
    return v2 * v1;
}
COL operator * (double v1, COL v2) {
    COL V = { static_cast<COL_t>(v2[0] * v1), static_cast<COL_t>(v2[1] * v1), static_cast<COL_t>(v2[2] * v1) };
    return V;
}
COL operator + (COL c1, COL c2) {
    return { static_cast<COL_t>(c1[0] + c2[0]), static_cast<COL_t>(c1[1] + c2[1]), static_cast<COL_t>(c1[2] + c2[2]) };
}
double abs(VEC &v1) {
    double V = v1 * v1;
    return std::sqrt(V);
}

VEC cross(VEC v1, VEC v2) {
    VEC V = { v1[1] * v2[2] - v1[2] * v2[1], v2[0] * v1[2] - v1[0] * v2[2], v1[0] * v2[1] - v2[0] * v1[1] };
    return V;
}

VEC normalize(VEC& v) {
    return (1 / abs(v)) * v;
}
VEC normalize(VEC&& v) {
    VEC v1 = (1 / abs(v)) * v;
    return v1;
}
/// <summary>
///  Rotation around x axe
/// </summary>
/// <param name="angle"></param>
/// <param name="v"></param>
/// <returns></returns>
VEC xRotate(double angle, VEC v) {
    MATR xRotation = { std::array<double, 3>{1, 0, 0},
        std::array<double,3>{0, std::cos(angle), -std::sin(angle)},
        std::array<double, 3>{ 0, std::sin(angle), std::cos(angle)}};
    return xRotation * v;
}
VEC yRotate(double angle, VEC v) {
    MATR yRotation = { std::array<double, 3>{std::cos(angle), 0, std::sin(angle)},
        std::array<double,3>{0, 1, 0},
        std::array<double, 3>{ -std::sin(angle), 0, std::cos(angle)}};
    return yRotation * v;
}
VEC xyRotate(double angle_x, double angle_y, VEC v) {
    MATR xRotation = { std::array<double, 3>{1, 0, 0},
        std::array<double,3>{0, std::cos(angle_x), -std::sin(angle_x)},
        std::array<double, 3>{ 0, std::sin(angle_x), std::cos(angle_x)} };
    MATR yRotation = { std::array<double, 3>{std::cos(angle_y), 0, std::sin(angle_y)},
        std::array<double,3>{0, 1, 0},
        std::array<double, 3>{ -std::sin(angle_y), 0, std::cos(angle_y)} };
    return xRotation * yRotation * v;
}

/// <summary>
/// Class for collecting all pixels data, maybe useless
/// </summary>
class Picture {
public:
    COL* arr;
    void set_color(int& i, int& j, COL& c) {
        arr[i + j * n1] = c;
    }
    COL get_color(int& i, int& j) {
        return arr[i + j * n1];
    }
    Picture() { arr = new COL[n1 * n2]; }
    ~Picture() { delete[] arr; }
};

/// <summary>
/// This class for discribe light
/// RAII done, but Inheritance...
/// </summary>
class LightObj {
protected:
    VEC* d;           // direction (for direct) or position (for point)
    char* type;
    double* intensity;
public:
    /// <summary>
    /// basic constructors (brand new data)
    /// </summary>
    /// <param name="v"></param>
    /// <param name="s"></param>
    /// <param name="i"></param>
    LightObj(VEC const& v, char const& s, double const& i) {
        d = new VEC; *d = v;
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
    VEC get_d() {
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
    COL* color;
    double* reflective; // how it matt or glossy (-1 - matt, 0-10000 - glossy)
    double* specular;   // reflections on/off and how effective
public:
    GenericObject(COL& c, double& refl, double& spec) {
        color = new COL;         *color = c;
        reflective = new double; *reflective = refl;
        specular = new double;   *specular = spec;
    };
    GenericObject(GenericObject& t) : GenericObject(*t.color, *t.reflective, *t.specular) {}
    GenericObject(GenericObject&& t) {
        color = t.color; t.color = nullptr;
        reflective = t.reflective; t.reflective = nullptr;
        specular = t.specular; t.specular = nullptr;
    }
    ~GenericObject() {
        delete color;
        delete reflective;
        delete specular;
    }
    /// <summary>
    /// Function calculates closest intercection for ray (actually line) and object 
    /// </summary>
    /// <param name="O - camera center in affine space"></param> 
    /// <param name="V - vector for the Intercecting ray"></param> 
    /// <param name="t_min - minimum distance for intercection"></param>
    /// <param name="t_max - maximum distance for intercection"></param>
    /// <param name="intersection - function output variable"></param>
    virtual void Intercections(VEC& O, VEC& V, double& t_min, double& t_max, double& intersection) = 0;
    virtual VEC get_norm(VEC P) = 0;
    virtual double get_reflective() = 0;
    virtual double get_specular() = 0;
    virtual COL get_color() = 0;
};

class SphereObj : public GenericObject {
    VEC* center;
    double* radius;
public:
    SphereObj(VEC v, double r, COL c, double spec, double refl) : GenericObject(c, refl, spec) {
        radius = new double; *radius = r;
        center = new VEC; *center = v;
    };
    SphereObj(VEC v, double r, COL c, double spec) : SphereObj(v, r, c, spec, 0.0) {};
    SphereObj(VEC v, double r, COL c) : SphereObj(v, r, c, -1.0) {};

    SphereObj(SphereObj const& t) : SphereObj(*t.center, *t.radius, *t.color, *t.specular, *t.reflective) {};

    SphereObj(SphereObj&& t) : GenericObject(t) {// may not work as intended
        center = t.center; t.center = nullptr;
        radius = t.radius; t.radius = nullptr;
    }

    ~SphereObj() {
        delete center;
        delete radius;
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

    virtual void Intercections(VEC& O, VEC& V, double& t_min, double& t_max, double& closest_t) override {
        VEC C = *center;
        double r = *radius;
        VEC OC = O - C;
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
        return *GenericObject::color;
    }
    virtual double get_reflective() override {
        return *GenericObject::reflective;
    }
    virtual double get_specular() override {
        return *GenericObject::specular;
    }
    VEC get_center() {
        return *center;
    }
    virtual VEC get_norm(VEC P) override {
        VEC N = P - *center;
        return N;
    }
};

class PlaneObj : public GenericObject {
private:
    VEC *norm;
    VEC *param;
public:
    PlaneObj(VEC norm, VEC param, COL c, double spec, double refl) : GenericObject(c, refl, spec) {
        PlaneObj::norm = new VEC; *PlaneObj::norm = norm;
        PlaneObj::param = new VEC; *PlaneObj::param = param;
    }
    PlaneObj(PlaneObj const& t) : PlaneObj(*t.norm, *t.param, *t.color, *t.specular, *t.reflective) {};
    ~PlaneObj() {
        delete norm;
        delete param;
    }

    PlaneObj& operator = (PlaneObj const& s) {
        PlaneObj tmp(s);
        std::swap(this->norm, tmp.norm);
        std::swap(this->param, tmp.param);
        std::swap(this->color, tmp.color);
        std::swap(this->reflective, tmp.reflective);
        std::swap(this->specular, tmp.specular);
        return *this;
    }

    virtual void Intercections(VEC& O, VEC& V, double& t_min, double& t_max, double& closest_t) override {
        double dot = *norm * V;
        double t1 = positive_inf;
        closest_t = positive_inf;
        if (abs(dot) >= epsilon) {
            VEC W = O-*param;
            double fac = -(*norm * W)/dot;
            if (fac > 0) {
                t1 = fac;
            }
        }
        if ((t1 >= t_min) && (t1 <= t_max)) {
            closest_t = t1;
        }
    }

    virtual VEC get_norm(VEC P) override {
        return *norm;
    }
    virtual COL get_color() override {
        return *GenericObject::color;
    }
    virtual double get_reflective() override {
        return *GenericObject::reflective;
    }
    virtual double get_specular() override {
        return *GenericObject::specular;
    }
};

class TriangleObj : public GenericObject {
private:
    VEC* norm;
    VEC* point_1;  //FIXED
    VEC* point_2;
    VEC* point_3;
public:
    TriangleObj(VEC normal, VEC Point_1, VEC Point_2, VEC Point_3, COL c, double spec, double refl): GenericObject(c, refl, spec) {
        norm = new VEC; *norm = normal;
        point_1 = new VEC; *point_1 = Point_1;
        point_2 = new VEC; *point_2 = Point_2;
        point_3 = new VEC; *point_3 = Point_3;
    }
    TriangleObj(TriangleObj const& t) : TriangleObj(*t.norm, *t.point_1, *t.point_2, *t.point_3, *t.color, *t.specular, *t.reflective) {};
    ~TriangleObj() {
        delete norm;
        delete point_1;
        delete point_2;
        delete point_3;
    }

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

    virtual void Intercections(VEC& O, VEC& V, double& t_min, double& t_max, double& closest_t) override {
        VEC e1 = *point_2 - *point_1;
        VEC e2 = *point_3 - *point_1;
        VEC pvec = cross(V, e2);
        double det = (e1 * pvec);
        double t = positive_inf;
        closest_t = positive_inf;
        if (!(det < epsilon && det > -epsilon)) {
            double inv_det = 1 / det;
            VEC tvec = O - *point_1;
            double u = (tvec * pvec) * inv_det;
            if (!(u < 0 || u > 1)) {
                VEC qvec = cross(tvec, e1);
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

    virtual VEC get_norm(VEC P) override {
        return *norm;
    }
    virtual COL get_color() override {
        return *GenericObject::color;
    }
    virtual double get_reflective() override {
        return *GenericObject::reflective;
    }
    virtual double get_specular() override {
        return *GenericObject::specular;
    }
};

class Render final {
    std::vector<GenericObject*> spheres = 
    {
         new SphereObj({0, 0, 15}, 2, BLUE, 500, 0.2),
         new SphereObj({4, -2, 20}, 1.5, RED, 100, 0),
         new PlaneObj({0, 1, 0}, {0, -3, 10}, {100, 100, 100}, 0.1, 0.1)
    };

    //std::vector<GenericObject*> spheres;
    std::vector<LightObj> lights = {
        LightObj('a', 0.3),
        LightObj({1, 1, -2}, 'd', 0.7)
    
    };

public:
    ~Render() { for (int i = 0; i < spheres.size(); i++) { delete spheres[i]; } }
    Render() {};
    /*
    Render(int n, GenericObject *obj1, ...) {  //FIXIT
        for (int i = 0; i < n; i++) {
            GenericObject* obj = obj1 + i * sizeof(obj1);
            spheres.push_back(obj);
            //obj1 = nullptr;
        }
    }
    */
    VEC ReflectRay(VEC R, VEC N) {
        // for sphere
        VEC V = 2 * (R * N) * N - R;
        return V;
    }
    void ClosestIntersection(VEC O, VEC V, double t_min, double t_max, GenericObject*& closest_sphere, double& closest_t) {
        for (int j = 0; j < spheres.size(); j++) {
            auto Sphere = spheres[j];
            double t;
            Sphere->Intercections(O, V, t_min, t_max, t);  // DONE
            if (((t >= t_min) && (t <= t_max)) && (t < closest_t)) {
                closest_t = t;
                closest_sphere = (spheres[j]);
            }
        }
    }

    /// <summary>
    /// </summary>
    /// <param name="P - the intersection point of the ray and the object"></param>
    /// <param name="N - normal for the surface in the point P"></param>
    /// <param name="V - vector from the P to the source of ray (point O)"></param>
    /// <param name="s - specularity of the object (матовость)"></param>
    /// <returns></returns>
    double ComputeLight(VEC P, VEC N, VEC V, double s) {
        // intensity of the initial color, after computing light
        double i = 0.0;

        // lights - vector<LightObj>
        for (int j = 0; j < lights.size(); j++) {
            LightObj Light = lights[j];
            if (Light.get_type() == 'a') {
                i += Light.get_intensity();
            }
            else {
                VEC L; // vector to the light
                if (Light.get_type() == 'p') {
                    L = Light.get_d() - P;
                }

                if (Light.get_type() == 'd') {
                    L = Light.get_d();
                }

                double closest_t = positive_inf;  // parameter of the distamce to the shadow_sphere
                GenericObject* shadow_sphere = nullptr;  // shadow_sphere - object, which will make a shadow
                ClosestIntersection(P, L, epsilon, positive_inf, shadow_sphere, closest_t); // calculate intersection
                
                if (shadow_sphere != nullptr) { // if the shadow_sphere != nullptr, then there is a shadow
                    continue;
                }
                // if not, then compute light
                double n_dot_l = N * L;
                if (n_dot_l > 0) {
                    if ((N*N != 0) && (L*L != 0)) {
                        i += Light.get_intensity() * n_dot_l / abs(N) / abs(L);
                    }
                }

                if (s != -1) {
                    VEC R = ReflectRay(L, N);
                    double r_dot_v = (R * V);
                    if ((r_dot_v > 0) && (R*R != 0) && (V*V != 0)) {
                        i += Light.get_intensity() * std::pow(r_dot_v / (abs(R) * abs(V)), s);
                    }
                }
            }

        }
        return i < 1 ? i : 1;
    }





    COL Trace(VEC O, VEC V, double t_min, double t_max, int depth) {
        double closest_t = positive_inf;
        GenericObject* closest_sphere = nullptr;

        // spheres (all of them)
        ClosestIntersection(O, V, t_min, t_max, closest_sphere, closest_t);

        if (closest_sphere == nullptr) {
            return BG_C;
        }

        VEC P = O - (-closest_t) * V;
        VEC N = closest_sphere->get_norm(P);

        if (abs(N) != 0) { N = (1 / abs(N)) * N; }

        COL local_color = ComputeLight(P, N, -V, closest_sphere->get_specular()) * closest_sphere->get_color();

        double r = closest_sphere->get_reflective();

        if ((depth <= 0) || (r <= 0)) {
            return local_color;
        }
        VEC R = ReflectRay(-V, N);
        COL reflected_color = Trace(P, R, epsilon, positive_inf, depth - 1);
        return (1 - r) * local_color + r * reflected_color;
    }
};
/*
int main() {
    Render r;
    Picture p;
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            VEC D = { (i - (double)n1 / 2) * (double)w / (double)n1,
                        -(j - (double)n2 / 2) * (double)h / (double)n2, (double)d };
            VEC O = { 0.0, 0.0, 0.0 };
            p.set_color(i, j, r.Trace(O, D, 1, positive_inf, 4));
        }
    }
    return 0;
}
*/