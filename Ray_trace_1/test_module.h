#include "framework.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
//#include<amp.h>
//#include <omp.h>

#ifdef OMP
int const MAX_THREADS = 3;
#endif

int const WIDTH = 30;
int const HEIGHT = 30;

int n1 = 300;  //  rename?
int n2 = 200;
int pixel = n1 / WIDTH;
double  d = 1;  // params of FOV and scaling screen // how does it work actually??
double  w = d;  // params of FOV and scaling screen
double h = (double)n2/n1 * d;  // params of FOV and scaling screen
double const positive_inf = 10000000;
double epsilon = 0.00001;
using VEC = std::array<double, 3>;
using COL_t = BYTE;
using COL = std::array<COL_t, 3>;
using MATR = std::array<VEC, 3>;

clock_t s_1 = 0;
clock_t e_1 = 0;

/// <summary>
/// Gram moment
/// </summary>
MATR Gram = 
{
    std::array<double, 3>{1, 0, 0},
    std::array<double, 3>{0, 1, 0},
    std::array<double, 3>{0, 0, 1}
};

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
#ifdef OMP
double operator * (VEC v1, VEC v2) {
    double answer = 0;
#pragma omp parallel for reduction(+:answer)
        for (int i = 0; i < 3; i++) {
            answer += v1[i] * v2[i];
        }
#pragma omp barrier
    return answer;
}
#else
double operator * (VEC v1, VEC v2) {
    //v2 = Gram * v2;  // epic boost
    double answer = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    return answer;
}
#endif

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
VEC xRotate(double angle, VEC v) {   // FIXIT
    MATR xRotation = { 
        VEC{1, 0, 0},
        VEC{0, std::cos(angle), -std::sin(angle)},
        VEC{ 0, std::sin(angle), std::cos(angle)}};
    return xRotation * v;
}
VEC yRotate(double angle, VEC v) {    // FIXIT
    MATR yRotation = { 
        VEC{std::cos(angle), 0, std::sin(angle)},
        VEC{0, 1, 0},
        VEC{ -std::sin(angle), 0, std::cos(angle)}};
    return yRotation * v;
}

VEC zRotate(double angle, VEC v) {    // FIXIT
    MATR zRotation = {
        VEC{std::cos(angle), -std::sin(angle), 0},
        VEC{std::sin(angle), std::cos(angle), 0},
        VEC{0, 0, 1} };
    return zRotation * v;
}
VEC xyRotate(double angle_x, double angle_y, VEC v) { // FIXIT
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
    /// <param name="R"></param>
    /// <param name="N"></param>
    /// <returns></returns>
VEC ReflectRay(VEC& R, VEC& N) {
    return  2 * (R * N) * N - R;
}

VEC ReflectRay(VEC&& R, VEC& N) {
    return  2 * (R * N) * N - R;
}

void PaintRect(HDC hdc, RECT* rect, COLORREF colour)
{
    COLORREF oldcr = SetBkColor(hdc, colour);
    ExtTextOut(hdc, 0, 0, ETO_OPAQUE, rect, (LPCWSTR)"", 0, 0);
    SetBkColor(hdc, oldcr);
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
    std::vector<GenericObject*>spheres;
    std::vector<LightObj> lights = {
        LightObj('a', 0.3),
        LightObj({1, 1, -2}, 'd', 0.3)
    };
public:
    // render camera motion parameters
    double Pi = 3.1415926536;
    VEC O = { 0, 0, 0 };                 // Ёто в общем точка, из которой лучи испускаютс€
    VEC DIR = normalize({ 0.1, 0, 1 });  // Looking direction
    VEC UP = normalize({ 0, 1, 0 });
    double step = 0.5;                   // how much you will move // for debug, all logic must be rewrite
    double angle = 0.005;                // how much you rotate

    // winAPI variables
    HINSTANCE hInst;                                // текущий экземпл€р
    WCHAR szTitle[MAX_LOADSTRING];                  // “екст строки заголовка
    WCHAR szWindowClass[MAX_LOADSTRING];            // им€ класса главного окна

    //~Render() { for (auto i : spheres) { delete i; } } // FIXIT
    Render(HINSTANCE hInstance) {
        hInst = hInstance;
        LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
        LoadStringW(hInstance, IDC_RAYTRACE1, szWindowClass, MAX_LOADSTRING);
    };

    Render operator = (Render r) {
        spheres.clear();
        for (auto i : r.spheres) {
            spheres.push_back(i);
        }
        return *this;
    }
    /// <summary>
    /// Use this function to set objects, that will be in render
    /// If you send object to render, it will now manage memory
    /// </summary>
    /// <param name="id"></param>
    /// <returns></returns>
    /// 
    void set_obj(GenericObject* object) {
        spheres.push_back(object);
    }
    
    GenericObject* get_object(int id) { // FIXIT (pointless)
        if (spheres.size() > id) {
            return spheres[id];
        }
        else {
            return nullptr;
        }
    }
   

    ///
    /// Register window class
    /// 
    ATOM MyRegisterClass(LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM))
    {
        WNDCLASSEXW wcex;

        wcex.cbSize = sizeof(WNDCLASSEX);

        wcex.style = CS_HREDRAW | CS_VREDRAW;
        wcex.lpfnWndProc = WndProc;
        wcex.cbClsExtra = 0;
        wcex.cbWndExtra = 0;
        wcex.hInstance = hInst;
        wcex.hIcon = LoadIcon(hInst, MAKEINTRESOURCE(IDI_RAYTRACE1));
        wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
        wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
        wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_RAYTRACE1);
        wcex.lpszClassName = szWindowClass;
        wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

        return RegisterClassExW(&wcex);
    }

//
//   ‘”Ќ ÷»я: InitInstance(HINSTANCE, int)
//
//   ÷≈Ћ№: —охран€ет маркер экземпл€ра и создает главное окно
//
//    ќћћ≈Ќ“ј–»»:
//
//        ¬ этой функции маркер экземпл€ра сохран€етс€ в глобальной переменной, а также
//        создаетс€ и выводитс€ главное окно программы.
//
    BOOL InitInstance(int nCmdShow)
    {
        HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
            CW_USEDEFAULT, 0, n1, n2, nullptr, nullptr, hInst, nullptr);
        ShowWindow(hWnd, nCmdShow);
        UpdateWindow(hWnd);
        return TRUE;
    }

    void DrawScene(HWND hWnd) {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);

        RECT rect;
        GetClientRect(hWnd, &rect);

        pixel = rect.right / WIDTH;
        n1 = rect.right - rect.right % pixel + pixel;
        n2 = rect.bottom - rect.bottom % pixel + pixel;
        double h = (double)n2 / n1 * d;


        HDC hmdc = CreateCompatibleDC(hdc);
        HBITMAP bit = CreateCompatibleBitmap(hdc, n1, n2);
        SelectObject(hmdc, bit);

        /*
        BITMAPINFO bif;
        ZeroMemory(&bif, sizeof(BITMAPINFO));
        bif.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bif.bmiHeader.biBitCount = 32;
        bif.bmiHeader.biWidth = n1;
        bif.bmiHeader.biHeight = n2;
        bif.bmiHeader.biPlanes = 1;
        bif.bmiHeader.biCompression = BI_RGB;
        bif.bmiHeader.biClrImportant = 0;
        */
#ifdef clock_1
        clock_t start = clock();
#endif
        COL c; VEC D;
        int i; int j;
        RECT l;

        for (int k = 0; k < n1 * n2; k += pixel) {
            i = k / n2 * pixel;
            j = k % n2;
            D = -cross(DIR, UP) * ((i - (double)n1 / 2) * w / (double)n1) - UP * ((j - (double)n2 / 2) * h / (double)n2) + DIR * d;
            c = Trace(O, D, 1, positive_inf, 0);  // set_color
            //HBRUSH hb = CreateSolidBrush(RGB(c[0], c[1], c[2]));
            l.left = i; l.top = j; l.right = i + pixel;  l.bottom = j + pixel;
            PaintRect(hmdc, &l, RGB(c[0], c[1], c[2]));
            //FillRect(hmdc, &l, hb); // this is huge
            //DeleteObject(hb);
            // SetPixel(hmdc, i, j, RGB(c[0], c[1], c[2]));
        }

        BitBlt(hdc, 0, 0, n1, n2, hmdc, 0, 0, SRCCOPY); // fast, 0ms

        //GetDIBits(hmdc, bit, 0, 0, 0, &bif, DIB_RGB_COLORS);
        //GetDIBits(hmdc, bit, 0, n2, im, &bif, DIB_RGB_COLORS);
        //SetDIBitsToDevice(hdc, 0, 0, n1, n2, 0, 0, 0, n2, im, &bif, DIB_RGB_COLORS);
        //SelectObject(hmdc, bit);

        //BitBlt(hdc, 0, 0, n1, n2, hmdc, 0, 0, SRCCOPY);
#ifdef clock_1
        clock_t end = clock();
        clock_t result = end - start;
#endif
        DeleteObject(bit);
        DeleteDC(hmdc);

        EndPaint(hWnd, &ps);
        UpdateWindow(hWnd);
    }

    void CameraRotate(HWND &hWnd) {
        POINT p;
        RECT rectangle;
        LPPOINT point = &p;
        LPRECT rect = &rectangle;
        GetCursorPos(point);
        GetWindowRect(hWnd, rect);
        int delta_y = point->x - (rect->right + rect->left) / 2;
        int delta_x = point->y - (rect->bottom + rect->top) / 2;
        if (abs(DIR[0]) < abs(DIR[2])) {
            if (DIR[2] > 0) {
                DIR = normalize(xRotate(delta_x * angle, DIR));
            }
            else {
                DIR = normalize(xRotate(-delta_x * angle, DIR));
            }
        }
        else {
            if (DIR[0] > 0) {
                DIR = normalize(zRotate(-delta_x * angle, DIR));
            }
            else {
                DIR = normalize(zRotate(delta_x * angle, DIR));
            }
        }
        DIR = normalize(yRotate(delta_y * angle, DIR));


        SetCursorPos((rect->right + rect->left) / 2, (rect->bottom + rect->top) / 2);
    }

    void ClosestIntersection(VEC O, VEC V, double t_min, double t_max, GenericObject*& closest_sphere, double& closest_t) {
        for (int j = 0; j < spheres.size(); j++) {
            auto Sphere = spheres[j];
            double t;
            Sphere->Intercections(O, V, t_min, t_max, t);  
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
                VEC L;                                                                     // vector to the light
                if (Light.get_type() == 'p') {                                             // FIXIT (refactor Light to be inheritance
                    L = Light.get_d() - P;
                }
                if (Light.get_type() == 'd') {
                    L = Light.get_d();
                }

                double closest_t = positive_inf;                                           // parameter of the distamce to the shadow_sphere
                GenericObject* shadow_sphere = nullptr;                                    // shadow_sphere - object, which will make a shadow
                ClosestIntersection(P, L, epsilon, positive_inf, shadow_sphere, closest_t);// calculate intersection
                
                if (shadow_sphere != nullptr) {                                            // if the shadow_sphere != nullptr, then there is a shadow
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
                        i += Light.get_intensity() * std::pow(r_dot_v / (abs(R) * abs(V)), s);// this is expensive FIXIT
                    }
                }
            }

        }
        return i < 1 ? i : 1;
    }

    COL Trace(VEC O, VEC V, double t_min, double t_max, int depth) {
        double closest_t = positive_inf;
        GenericObject* closest_sphere = nullptr;

        if (spheres.size() > 0) 
        {
            ClosestIntersection(O, V, t_min, t_max, closest_sphere, closest_t);
        }

        if (closest_sphere == nullptr) 
        {
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

        COL reflected_color = Trace(P, ReflectRay(-V, N), epsilon, positive_inf, depth - 1);
        return (1 - r) * local_color + r * reflected_color;
    }
};
