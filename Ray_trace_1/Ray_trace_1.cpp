﻿// Ray_trace_1.cpp : Определяет точку входа для приложения.

#define _SILENCE_AMP_DEPRECATION_WARNINGS
#define debug
#define MAX_LOADSTRING 100
#define PICTURE 10

#include "Ray_trace_1.h"
#include "test_module.h"
#include <time.h>

LRESULT CALLBACK  WndProc(HWND, UINT, WPARAM, LPARAM);

Render r; // Рендер ( это кринж )

double box_param = 20;
double coridor_param = 2;

#ifdef GPU
std::vector<COL_t> array_pixel(n1*n2*3);
concurrency::array<COL_t, 2> array_pixel_gpu(n1* n2, 3, array_pixel.begin(), array_pixel.end());
//concurrency::array_view<COL_t, 2> array_pixel_gpu(n1 * n2, 3, array_pixel);
#endif 

#if PICTURE == 10
// TODO : realization
void draw_planets() {
    
}
#endif

/// <summary>
/// main function of the program
/// </summary>
/// <param name="hInstance"></param>
/// <param name="hPrevInstance"></param>
/// <param name="lpCmdLine"></param>
/// <param name="nCmdShow"></param>
/// <returns></returns>
int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
    _In_opt_ HINSTANCE hPrevInstance,
    _In_ LPWSTR    lpCmdLine,
    _In_ int       nCmdShow)
{
    // init
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    r = Render(hInstance); // Рендер (инициализировали)
#if PICTURE == 404
    depth = 5;
    r.set_position({ -8, 2, -8 });
    r.set_direction({ 1, 0, 4});
    r.set_obj(new SphereObj({ 0.7, 1, 13 }, 1.0, GREENYELLOW, 5, 0.2));
    r.set_obj(new SphereObj({ -0.7, 1, 13 }, 1.0, GREENYELLOW, 5, 0.2));
    r.set_obj(new SphereObj({ 0, 2.0, 13 }, 1.0, GREENYELLOW, 5, 0.2));
    r.set_obj(new SphereObj({ 0, 3.0, 13 }, 1.0, GREENYELLOW, 5, 0.2));
    r.set_obj(new SphereObj({ 0, 4.0, 13 }, 1.1, RED, 100, 0.02));
    r.set_obj(new WallObj({ box_param, 0, box_param}, { box_param, 10, -box_param }, LIGHTCORAL, 10, 0.01));
    r.set_obj(new WallObj({ box_param, 0, box_param }, { -box_param, 10, box_param }, LIGHTCORAL, 10, 0.01));
    r.set_obj(new WallObj({-box_param, 0, -box_param }, { -box_param, 10, box_param }, LIGHTCORAL, 10, 0.01));
    r.set_obj(new WallObj({-box_param, 0, -box_param }, { box_param, 10, -box_param }, LIGHTCORAL, 10, 0.01));
    r.set_obj(new WallObj({-9, 1, box_param-0.01 },{9, 9, box_param-0.01 }, {2, 2, 2}, 10, 0.95));
    r.set_obj(new PlaneObj({ 0, -1, 0 }, { 0, 0, 0}, { 200, 210, 210 }, 10, 0));
    r.set_light(new LightObj('a', 0.5));
    r.set_light(new LightObj({ 0, 5, -1 }, 'd', 0.3));
    
#endif
#if PICTURE == 1
    BG_C = { 10, 30, 140 };
    depth = 5;
    r.set_position({ -40, 2, 0 });
    r.set_direction({ 3, 0, 1 });
    r.set_obj(new PlaneObj({ 0, -1, 0 }, { 0, -1, 0 }, SADBROWN, 100, 0.01));
    r.set_obj(new PlaneObj({ 0, 0, -1 }, { 0, 0, 20 }, { 0, 0, 0 }, 100, 0.95));
    r.set_light(new LightObj({ -10, 10, -10 }, 'd', 0.4));
    r.set_light(new LightObj('a', 0.6));
    r.set_obj(new SphereObj({ 0, 5.0, 13 }, 4.0, LIGHTCORAL, -1, 0.1));
#endif
#if PICTURE == 10
    BG_C = { 0, 0, 0 };
    depth = 5;
    r.set_position({ 0, 4, -20 });
    r.set_direction({ 0, -1, 4 });
    r.set_obj(new SphereObj({ 0, 0, 0 }, 4.0, { 200, 200, 0 }, -1, 0));
    r.set_obj(new SphereObj({ 6, 0, 0 }, 0.8, { 0, 200, 0 }, 100, 0.1));
    r.set_obj(new SphereObj({ -7, 0, -1 }, 1.5, LIGHTCORAL, 100, 0.5));
    r.set_obj(new SphereObj({ -4, 0.4, -9 }, 0.5, { 200, 180, 200 }, 10, 0.1));
    r.set_obj(new SphereObj({ 8, 0.8, 10 }, 0.5, { 0, 180, 200 }, 10, 0.1));
    r.set_light(new LightObj('a', 0.6));
    r.set_light(new LightObj({ -10, 10, -10 }, 'd', 0.4));
#endif
#ifdef OMP
    omp_set_num_threads(MAX_THREADS);
#endif

    // Выполнить инициализацию приложения:
    r.MyRegisterClass(&WndProc); // так нужно, передаём ссылку на функцию, чего таково
    r.InitInstance(nCmdShow);

    // переменная для сообщений, и переменная для (чего?)
    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_RAYTRACE1));
    MSG msg;

    // Цикл основного сообщения 
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg); // диспетчеризация событий
        }
    }
    r.delete_spheres();
    return (int) msg.wParam;
}

// I think this will never happen
#ifdef GPU
void DrawOptimal(HWND hWnd) {
    PAINTSTRUCT ps;
    HDC hdc = BeginPaint(hWnd, &ps);

    RECT rect;
    GetClientRect(hWnd, &rect);

    HDC hmdc = CreateCompatibleDC(hdc);
    HBITMAP bit = CreateCompatibleBitmap(hdc, n1, n2);
    SelectObject(hmdc, bit);


    concurrency::parallel_for_each(
        array_pixel_gpu.extent, [=](concurrency::index<2> idx) restrict(amp) {

            //auto i = idx[0] / n2;
            //auto j = idx[0] % n2;
            //VEC D = -cross(DIR, UP) * ((i - (double)n1 / 2) * w / (double)n1) - UP * ((j - (double)n2 / 2) * h / (double)n2) + DIR * d;

        }
    );
    array_pixel = array_pixel_gpu;
#ifdef clock_1
    clock_t start = clock();
#endif
    // Magic happens here


    BitBlt(hdc, 0, 0, n1, n2, hmdc, 0, 0, SRCCOPY); // fast, 0ms
#ifdef clock_1
    clock_t end = clock();
#endif
    DeleteObject(bit);
    DeleteDC(hmdc);

    EndPaint(hWnd, &ps);
    UpdateWindow(hWnd);
}
#endif

//
//  ФУНКЦИЯ: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  ЦЕЛЬ: Обрабатывает сообщения в главном окне.
//
//  WM_COMMAND  - обработать меню приложения
//  WM_PAINT    - Отрисовка главного окна
//  WM_DESTROY  - отправить сообщение о выходе и вернуться
//
// dispatcherisation function
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Разобрать выбор в меню:
            switch (wmId)
            {
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_KEYDOWN:
    {
        WPARAM key = wParam; //Получаем код нажатой клавиши

        // x translation
        if (key == 65) {
            r.O = r.O - r.step * normalize(cross(UP, r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 68) {
            r.O = r.O + r.step * normalize(cross(UP, r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }

        // y translation
        if (key == 83) {
            r.O = r.O + r.step * normalize(cross(cross(UP, r.DIR), r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 87) {
            r.O = r.O - r.step * normalize(cross(cross(UP, r.DIR), r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }

        // z translation
        if (key == 81) {
            r.O = r.O + r.step * r.DIR;
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 69) {
            r.O = r.O - r.step * r.DIR;
            InvalidateRect(hWnd, NULL, NULL);
        }

        // DIR change
        if (key == 37) {
            r.DIR = normalize(yRotate(-r.angle, r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 39) {
            r.DIR = normalize(yRotate(r.angle, r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 38) {
            r.DIR = normalize(xRotate(-r.angle, r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 40) {
            r.DIR = normalize(xRotate(r.angle, r.DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }

    }
    /*
    case WM_MOUSEMOVE:
    {
        r.CameraRotate(hWnd);

        InvalidateRect(hWnd, NULL, NULL);
    }
    break;
    */
    case WM_PAINT:
    {
        r.DrawScene(hWnd);
    }
    break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}
