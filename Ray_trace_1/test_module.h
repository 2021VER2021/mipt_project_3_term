#include <iostream>
#include <vector>
#include <cmath>
#include <array>

int const n1 = 800;  //  rename?
int const n2 = 400;
double const d = 1;  // params of FOV and scaling screen // how does it work actually??
double const w = d;  // params of FOV and scaling screen
double const h = (double)n2/n1 * d;  // params of FOV and scaling screen
double const positive_inf = 100000000;
double epsilon = 0.00001;
using VEC = std::array<double, 3>;
using COL_t = unsigned char;
using COL = std::array<COL_t, 3>;
using MATR = std::array<VEC, 3>;

VEC v1 = { 1, 0, 0 };
VEC v2 = { 0, 1, 0 };
VEC v3 = { 0, 0, 1 };
MATR Gramm = { v1, v2, v3 };

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
double operator * (VEC v1, VEC v2) {
    double V = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            V += v1[i] * v2[j] * Gramm[i][j];
        }
    }
    return V;
}
VEC operator * (double v1, VEC v2) {
    VEC V = { v2[0] * v1, v2[1] * v1, v2[2] * v1 };
    return V;
}
COL operator * (double v1, COL v2) {
    COL V = { static_cast<COL_t>(v2[0] * v1), static_cast<COL_t>(v2[1] * v1), static_cast<COL_t>(v2[2] * v1) };
    return V;
}
COL operator + (COL c1, COL c2) {
    return { static_cast<COL_t>(c1[0] + c2[0]), static_cast<COL_t>(c1[1] + c2[1]), static_cast<COL_t>(c1[2] + c2[2]) };
}
double abs(VEC v1) {
    double V = v1 * v1;
    return std::sqrt(V);
}
/// <summary>
/// Class for collecting all pixels data, maybe useless
/// </summary>
class Picture {
public:
    COL* arr;
    void set_color(int i, int j, COL c) {
        arr[i + j * n1] = c;
    }
    COL get_color(int i, int j) {
        return arr[i + j * n1];
    }
    Picture() { arr = new COL[n1 * n2]; }
    ~Picture() { delete[] arr; }
};

/// <summary>
/// This class for discribe light
/// RAII done
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
    /// <param name="O - coord center in affine space"></param> 
    /// <param name="V - vector for the Intercecting ray"></param> 
    /// <param name="t_min - minimum distance for intercection"></param>
    /// <param name="t_max - maximum distance for intercection"></param>
    /// <param name="intersection - function output variable"></param>
    virtual void Intercections(VEC O, VEC V, double t_min, double t_max, double& intersection) = 0;
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

    virtual void Intercections(VEC O, VEC V, double t_min, double t_max, double& closest_t) override {
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

    virtual void Intercections(VEC O, VEC V, double t_min, double t_max, double& closest_t) override {
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

class Render final {
    std::vector<GenericObject*> spheres = { new SphereObj({-0.7, -0.5, 32}, 1, RED),
                     new SphereObj({0.7, -0.5, 32}, 1, RED),
                     new SphereObj({0, 0.6, 32}, 1, YELLOW, 10, 0.2),
                     new SphereObj({0, 1.7, 32}, 1, {116, 66, 200}, 500, 0.1),
        new SphereObj({0, 2.2, 32}, 1.1, BLUE, 500, 0.5),
    new PlaneObj({0, 1, 0}, {0, -3, 10}, {255, 255, 255}, -1, 0.7)};

    //std::vector<GenericObject*> spheres;
    std::vector<LightObj> lights = { LightObj('a', 0.3),
                                 LightObj({1, 1, -2}, 'd', 0.7) };

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
        // spheres
        for (int j = 0; j < spheres.size(); j++) {
            auto Sphere = spheres[j];
            double t;
            Sphere->Intercections(O, V, t_min, t_max, t);  // DONE
            if (((t >= t_min) && (t <= t_max)) && (t < closest_t)) {
                closest_t = t;
                closest_sphere = (spheres[j]);
            }
        }
        // more code here to intersect other objects? or something
    }
    double ComputeLightDop(VEC P, VEC N, VEC V, double s) {
        double i = 0.0;
        for (int j = 0; j < lights.size(); j++) {
            LightObj Light = lights[j];
            if (Light.get_type() == 'a') {
                i += Light.get_intensity();
            }
            else {
                VEC L;
                if (Light.get_type() == 'p') {
                    L = Light.get_d() - P;
                }
                else {
                    L = Light.get_d();
                }
                double closest_t = positive_inf;
                GenericObject* shadow_sphere = nullptr;

                // spheres 
                ClosestIntersection(P, L, epsilon, positive_inf, shadow_sphere, closest_t);
                if (shadow_sphere != nullptr) {
                    continue;
                }

                double n_dot_l = N * L;
                if (n_dot_l > 0) {
                    if ((abs(N) != 0) && (abs(L) != 0)) {
                        i += Light.get_intensity() * n_dot_l / abs(N) / abs(L);
                    }
                }

                if (s != -1) {
                    VEC R = ReflectRay(L, N);
                    double r_dot_v = (R * V);
                    if ((r_dot_v > 0) && (abs(R) != 0) && (abs(V) != 0)) {
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

        COL local_color = ComputeLightDop(P, N, -V, closest_sphere->get_specular()) * closest_sphere->get_color();

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