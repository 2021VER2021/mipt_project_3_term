#include <iostream>
#include <vector>
#include <cmath>
#include <array>

int const n1 = 800;  //  rename?
int const n2 = 400;
double const d = 1;  // params of FOV and scaling screen
double const w = d;  // params of FOV and scaling screen
double const h = 0.5 * d;  // params of FOV and scaling screen
double const positive_inf = 10000000;
using VEC = std::array<double, 3>;
using COL_t = unsigned char;
using COL = std::array<COL_t, 3>;

const COL RED = { 250, 0, 0 };    // Map? maybe
const COL BLUE = { 0, 0, 250 };
const COL YELLOW = { 250, 250, 0 };
const COL BG_C = { 0, 191, 255};

VEC operator - (VEC v1, VEC v2) {
    VEC V = { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
    return V;
}
VEC operator - (VEC v1) {
    VEC V = { -v1[0], -v1[1], -v1[2]};
    return V;
}
double operator * (VEC v1, VEC v2) {
    double V = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
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
    double V = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
    return std::sqrt(V);
}

class Picture {
private:
    COL* arr;
public:
    void set_color(int i, int j, COL c) {
        arr[i + j * n1] = c;
    }
    COL get_color(int i, int j) {
        return arr[i + j * n1];
    }
    Picture() { arr = new COL[n1 * n2]; }
    ~Picture() { delete [] arr; }
};

struct LightObj {
    VEC d;           // direction (for direct) or position (for point)
    char type;
    double intensity;

    LightObj(VEC const& v, char const&  s, double const& i) : d(v), type(s), intensity(i) {};
    LightObj(LightObj const & obj): LightObj(obj.d, obj.type, obj.intensity) {};
    LightObj(char s, double i) : LightObj({ 0, 0, 0 }, s, i) {};

    LightObj& operator = (LightObj const &s) {
        d = s.d;
        type = s.type;
        intensity = s.intensity;
        return *this;
    }
};

class GenericObject{
protected:
    COL const color;
    double const reflective; // how it matt or glossy (-1 - matt, 0-10000 - glossy)
    double const specular;   // reflections on/off and how effective
public:
    virtual void Intercections(VEC O, VEC V, double arr[]) = 0;
    GenericObject(COL color, double reflective, double specular) : color(color), reflective(reflective), specular(specular) {};
    virtual VEC get_norm(VEC P) = 0;
    virtual double get_reflective() = 0;
    virtual double get_specular() = 0;
    virtual COL get_color() = 0;
};

class SphereObj final: public GenericObject{
    VEC const center;
    double const radius;
public:
    SphereObj(VEC v, double r, COL c, double spec, double refl) : center(v), radius(r), GenericObject(c, refl, spec) {};
    SphereObj(VEC v, double r, COL c, double spec) : SphereObj(v, r, c, spec, 0) {};
    SphereObj(VEC v, double r, COL c) : SphereObj(v, r, c, -1) {};
    ~SphereObj(){
    
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

    VEC get_center() {
        return center;
    }

    virtual void Intercections(VEC O, VEC V, double intersection[]) override {
        VEC C = center;
        double r = radius;
        VEC OC = O - C;
        double k1 = V * V;
        double k2 = 2 * (OC * V);
        double k3 = (OC * OC) - r * r;
        double diskr = k2 * k2 - 4 * k1 * k3;
        if (diskr < 0) {
            intersection[0] = positive_inf;
            intersection[1] = positive_inf;
            return;
        }
        intersection[0] = (-k2 + std::sqrt(diskr)) / (2 * k1);
        intersection[1] = (-k2 - std::sqrt(diskr)) / (2 * k1);
    };

    virtual VEC get_norm(VEC P) override {
        VEC N = P - center;
        return N;
    }
};

class Render final{
    std::vector<GenericObject*> spheres = {new SphereObj({0, -0.2, 8}, 1, RED, -1, 0.95),
                     new SphereObj({2, 0, 5}, 1, BLUE),
                     new SphereObj({0, -5001, 0}, 5000, YELLOW, 10),
                     new SphereObj({-2, 0, 6}, 1, {116, 66, 200}, 500, 0.1) };
    
    //std::vector<GenericObject*> spheres;
    std::vector<LightObj> lights = { LightObj('a', 0.3),
                                 LightObj({1, 1, -2}, 'd', 0.7) };

public:
    ~Render() {for (int i = 0; i < spheres.size(); i++){delete spheres[i];} }
    Render(){};
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
            double arr[2];
            Sphere->Intercections(O, V, arr);  // DONE
            double t1 = arr[0];
            double t2 = arr[1];
            if (((t1 >= t_min) and (t1 <= t_max)) and (t1 < closest_t)) {
                closest_t = t1;
                closest_sphere = (spheres[j]);
            }
            if (((t2 >= t_min) and (t2 <= t_max)) and (t2 < closest_t)) {
                closest_t = t2;
                closest_sphere = (spheres[j]);
            }
        }
        // more code here to intersect other objects? or something
    }
    double ComputeLightDop(VEC P, VEC N, VEC V, double s) {
        double i = 0.0;
        for (int j = 0; j < lights.size(); j++) {
            LightObj Light = lights[j];
            if (Light.type == 'a') {
                i += Light.intensity;
            }
            else {
                VEC L;
                if (Light.type == 'p') {
                    L = Light.d - P;
                }
                else {
                    L = Light.d;
                }
                double closest_t = positive_inf;
                GenericObject* shadow_sphere = nullptr;

                // spheres 
                ClosestIntersection(P, L, 0.001, positive_inf, shadow_sphere, closest_t);
                if (shadow_sphere != nullptr) {
                    continue;
                }

                double n_dot_l = N * L;
                if (n_dot_l > 0) {
                    if ((abs(N) != 0) and (abs(L) != 0)) {
                        i += Light.intensity * n_dot_l / abs(N) / abs(L);
                    }
                }

                if (s != -1) {
                    VEC R = ReflectRay(L, N);
                    double r_dot_v = (R * V);
                    if ((r_dot_v > 0) and (abs(R) != 0) and (abs(V) != 0)) {
                        i += Light.intensity * std::pow(r_dot_v / (abs(R) * abs(V)), s);
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

        if ((depth <= 0) or (r <= 0)) {
            return local_color;
        }
        VEC R = ReflectRay(-V, N);
        COL reflected_color = Trace(P, R, 0.0001, positive_inf, depth - 1);
        return (1 - r) * local_color + r * reflected_color;
    }
};

/*
int main() {
    COL** pixels = new COL * [n1];
    for (int i = 0; i < n1; i++) {
        pixels[i] = new COL[n2];
    }

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            VEC D = { (i - (double)n1 / 2) * (double)w / n1, (j - (double)n2 / 2) * (double)h / n2, d };
            VEC O = { 0.0, 0.0, 0.0 };
            ///std::vector<double> P = V - O;
            pixels[i][j] = Trace(O, D, 1, positive_inf);
        }
    }


    for (int i = 0; i < n1; i++) {
        delete[] pixels[i];
    }
    delete[] pixels;
}
*/