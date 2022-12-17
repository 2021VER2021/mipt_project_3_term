#pragma once
#define MODEL_SEND 1024
#define MAX_LOADSTRING 100

#include "Resource.h"
#include "objects.h"

class Render final {
    std::vector<GenericObject*>objects;
    std::vector<LightObj> lights;
public:
    // render camera motion parameters
    double Pi = 3.1415926536;
    VEC3 O = { 0, 1, 0 };                 // Точка, из которой лучи испускаются
    VEC3 DIR = normalize({ 0, 0, 1 });  // Looking direction
    double step = 0.5;                   // how much you will move 
    double angle = 0.005;                // how much you rotate

    // winAPI variables
    HINSTANCE hInst;                                // текущий экземпляр
    WCHAR szTitle[MAX_LOADSTRING];                  // Текст строки заголовка
    WCHAR szWindowClass[MAX_LOADSTRING];            // имя класса главного окна

    //~Render() { for (auto i : objects) { delete i; } } // FIXIT
    Render() : hInst(nullptr) {};
    Render(HINSTANCE hInstance) {
        hInst = hInstance;
        LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
        LoadStringW(hInstance, IDC_RAYTRACE1, szWindowClass, MAX_LOADSTRING);
    };

    /// <summary>
    /// Use this function to set objects, that will be in render
    /// If you send object to render, it will now manage memory
    /// </summary>
    /// <param name="id"></param>
    /// <returns></returns>
    /// 
    void set_obj(GenericObject* object) {
        objects.push_back(object);
    }

    void set_light(LightObj* object) {
        lights.push_back(*object);
    }

    void set_position(VEC3 v) {
        O = v;
    }

    void set_direction(VEC3 v) {
        DIR = normalize(v);
    }

    GenericObject* get_object(int id) {
        if (objects.size() > id) {
            return objects[id];
        }
        else {
            return nullptr;
        }
    }

    void delete_objects() {
        for (auto i : objects) {
            delete i;
        }
        objects.clear();
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
    //   ФУНКЦИЯ: InitInstance(HINSTANCE, int)
    //
    //   ЦЕЛЬ: Сохраняет маркер экземпляра и создает главное окно
    //
    //   КОММЕНТАРИИ:
    //
    //        В этой функции маркер экземпляра сохраняется в глобальной переменной, а также
    //        создается и выводится главное окно программы.
    //
    BOOL InitInstance(int nCmdShow)
    {
        HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
            CW_USEDEFAULT, 0, n1, n2, nullptr, nullptr, hInst, nullptr);
        ShowWindow(hWnd, nCmdShow);
        UpdateWindow(hWnd);
        PostMessage(hWnd, MODEL_SEND, 0, 0);
        return TRUE;
    }

    void DrawScene(HWND hWnd) {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);

        RECT rect;
        GetClientRect(hWnd, &rect);

        //pixel = rect.right / WIDTH;
        //n1 = rect.right - rect.right % pixel + pixel;
        //n2 = rect.bottom - rect.bottom % pixel + pixel;
        //double h = (double)n2 / n1 * d;


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
        COL c; VEC3 D;
        int i; int j;
        RECT l;

        for (int k = 0; k < n1 * n2; k += pixel) {
            i = k / n2 * pixel;
            j = k % n2;
            D = -cross(DIR, UP) * ((i - (double)n1 / 2) * w / (double)n1) - UP * ((j - (double)n2 / 2) * h / (double)n2) + DIR * d;
            c = Trace(O, D, 0.01, positive_inf, depth);  // set_color
            //HBRUSH hb = CreateSolidBrush(RGB(c[0], c[1], c[2]));
            l.left = i; l.top = j; l.right = i + pixel;  l.bottom = j + pixel;
            PaintRect(hmdc, &l, RGB(c[0], c[1], c[2]));
            //FillRect(hmdc, &l, hb); // this is huge
            //DeleteObject(hb);
            //SetPixel(hmdc, i, j, RGB(c[0], c[1], c[2]));
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

    void CameraRotate(HWND& hWnd)
    {
        POINT p;
        RECT rectangle;
        LPPOINT point = &p;
        LPRECT rect = &rectangle;
        GetCursorPos(point);
        GetWindowRect(hWnd, rect);
        int delta_y = point->x - (rect->right + rect->left) / 2;
        int delta_x = point->y - (rect->bottom + rect->top) / 2;
        if (abs(DIR.x) < abs(DIR.z)) {
            if (DIR.z > 0) {
                DIR = normalize(xRotate(delta_x * angle, DIR));
            }
            else {
                DIR = normalize(xRotate(-delta_x * angle, DIR));
            }
        }
        else {
            if (DIR.x > 0) {
                DIR = normalize(zRotate(-delta_x * angle, DIR));
            }
            else {
                DIR = normalize(zRotate(delta_x * angle, DIR));
            }
        }
        DIR = normalize(yRotate(delta_y * angle, DIR));


        SetCursorPos((rect->right + rect->left) / 2, (rect->bottom + rect->top) / 2);
    }

    void ClosestIntersection(VEC3 O, VEC3 V, double t_min, double t_max, GenericObject*& closest_object, double& closest_t) {
        double t;
        GenericObject* object;
        for (int j = 0; j < objects.size(); j++) {
            object = objects[j];
            object->Intercections(O, V, t_min, t_max, t);
            if (((t >= t_min) && (t <= t_max)) && (t < closest_t)) {
                closest_t = t;
                closest_object = (objects[j]);
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
    double ComputeLight(VEC3& P, VEC3& N, VEC3 V, double s) {
        // intensity of the initial color, after computing light
        double i = 0.0;

        // lights - vector<LightObj>
        for (LightObj& Light : lights) {
            if (Light.get_type() == 'a') {
                i += Light.get_intensity();
            }
            else {
                VEC3 L;                                                                     // vector to the light
                if (Light.get_type() == 'p') {                                             // FIXIT (refactor Light to be inheritance
                    L = Light.get_d() - P;
                }
                if (Light.get_type() == 'd') {
                    L = Light.get_d();
                }

                double closest_t = positive_inf;                                           // parameter of the distamce to the shadow_object
                GenericObject* shadow_object = nullptr;                                    // shadow_object - object, which will make a shadow
                ClosestIntersection(P, L, ray_epsilon, positive_inf, shadow_object, closest_t);// calculate intersection

                if (shadow_object != nullptr) {                                            // if the shadow_object != nullptr, then there is a shadow
                    continue;
                }
                // if not, then compute light
                double n_dot_l = N * L;
                if (n_dot_l > 0) {
                    if ((N * N != 0) && (L * L != 0)) {
                        i += Light.get_intensity() * n_dot_l / abs(N) / abs(L);
                    }
                }

                if (s != -1) {
                    VEC3 R = ReflectRay(L, N);
                    double r_dot_v = (R * V);
                    if ((r_dot_v > 0) && (R * R != 0) && (V * V != 0)) {
                        i += Light.get_intensity() * std::pow(r_dot_v / (abs(R) * abs(V)), s);// this is expensive FIXIT
                    }
                }
            }

        }
        return i < 1 ? i : 1;
    }

    COL Trace(VEC3 O, VEC3 V, double t_min, double t_max, int depth) {
        double closest_t = positive_inf;
        GenericObject* closest_object = nullptr;

        if (objects.size() > 0)
        {
            ClosestIntersection(O, V, t_min, t_max, closest_object, closest_t);
        }

        if (closest_object == nullptr)
        {
            return BG_C;
        }

        VEC3 P = O + closest_t * V;
        VEC3 N = closest_object->get_norm(P);

        if (abs(N) != 0) { N = normalize(N); }

        COL local_color = ComputeLight(P, N, -V, closest_object->get_specular()) * closest_object->get_color();

        double r = closest_object->get_reflective();

        if ((depth <= 0) || (r <= 0)) {
            return local_color;
        }

        COL reflected_color = Trace(P, ReflectRay(-V, N), epsilon, positive_inf, depth - 1);
        return (1 - r) * local_color + r * reflected_color;
    }
};
