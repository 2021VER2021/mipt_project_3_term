// Ray_trace_1.cpp : Определяет точку входа для приложения.
#define _SILENCE_AMP_DEPRECATION_WARNINGS
#include "framework.h"
#include "Ray_trace_1.h"
//#include "ray_trace.h"
#include "test_module.h"
//#include "kernel.cu"
#include <time.h>
//#include <omp.h>
//#include <boost/compute.hpp>
#include<amp.h>

#define debug
#define MAX_LOADSTRING 100

// Глобальные переменные:
HINSTANCE hInst;                                // текущий экземпляр
WCHAR szTitle[MAX_LOADSTRING];                  // Текст строки заголовка
WCHAR szWindowClass[MAX_LOADSTRING];            // имя класса главного окна
Render r; // Рендер

VEC O = { 0, 0, 0 }; // Это в общем точка, из которой лучи испускаются
VEC DIR = normalize({ 0.1, 0, 1 });  // Looking direction
VEC UP = normalize({ 0, 1, 0 });
double step = 0.5; // how much you will move // for debug, all logic must be rewrite
double angle = 0.01; // how much you rotate
// Отправить объявления функций, включенных в этот модуль кода:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Разместите код здесь.

    // Инициализация глобальных строк
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_RAYTRACE1, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Выполнить инициализацию приложения:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_RAYTRACE1));

    MSG msg;

    // Цикл основного сообщения:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}

void Draw(HWND hWnd) {
    PAINTSTRUCT ps;
    HDC hdc = BeginPaint(hWnd, &ps);
    RECT rect;
    GetClientRect(hWnd, &rect);
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
    clock_t start = clock();
    COL c; VEC D;
    int i; int j;

    
    for (int k = 0; k < n1 * n2; k += pixel) {
        i = k / n2 * pixel;
        j = k % n2;
        D = -cross(DIR, UP) * ((i - (double)n1 / 2) * w / (double)n1) - UP * ((j - (double)n2 / 2) * h / (double)n2) + DIR * d;
        c = r.Trace(O, D, 1, positive_inf, 1);  // set_color
        HBRUSH hb = CreateSolidBrush(RGB(c[0], c[1], c[2]));
        RECT l;  l.left = i; l.top = j; l.right = i + pixel;  l.bottom = j + pixel;
        FillRect(hmdc, &l, hb);
        DeleteObject(hb);
        // SetPixel(hmdc, i, j, RGB(c[0], c[1], c[2]));
    }

    BitBlt(hdc, 0, 0, n1, n2, hmdc, 0, 0, SRCCOPY);

    //GetDIBits(hmdc, bit, 0, 0, 0, &bif, DIB_RGB_COLORS);
    //GetDIBits(hmdc, bit, 0, n2, im, &bif, DIB_RGB_COLORS);
    //SetDIBitsToDevice(hdc, 0, 0, n1, n2, 0, 0, 0, n2, im, &bif, DIB_RGB_COLORS);
    //SelectObject(hmdc, bit);

    //BitBlt(hdc, 0, 0, n1, n2, hmdc, 0, 0, SRCCOPY);

    clock_t end = clock();
    DeleteObject(bit);
    DeleteDC(hmdc);

    EndPaint(hWnd, &ps);
    UpdateWindow(hWnd);
}

//
//  ФУНКЦИЯ: MyRegisterClass()
//
//  ЦЕЛЬ: Регистрирует класс окна.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_RAYTRACE1));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_RAYTRACE1);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

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
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // Сохранить маркер экземпляра в глобальной переменной

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, n1, n2, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

//
//  ФУНКЦИЯ: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  ЦЕЛЬ: Обрабатывает сообщения в главном окне.
//
//  WM_COMMAND  - обработать меню приложения
//  WM_PAINT    - Отрисовка главного окна
//  WM_DESTROY  - отправить сообщение о выходе и вернуться
//
//
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
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
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
            O = O - step * normalize(cross(UP, DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 68) {
            O = O + step * normalize(cross(UP, DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        // y translation
        if (key == 83) {
            O = O - step * normalize(cross(cross(UP, DIR), DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 87) {
            O = O + step * normalize(cross(cross(UP, DIR), DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        // z translation
        if (key == 81) {
            O = O + step * DIR;
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 69) {
            O = O - step * DIR;
            InvalidateRect(hWnd, NULL, NULL);
        }
        // DIR change
        if (key == 37) {
            DIR = normalize(yRotate(-angle, DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 39) {
            DIR = normalize(yRotate(angle, DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 38) {
            DIR = normalize(xRotate(-angle, DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
        if (key == 40) {
            DIR = normalize(xRotate(angle, DIR));
            InvalidateRect(hWnd, NULL, NULL);
        }
    }
    case WM_MOUSEMOVE:
    {
        POINT p;
        RECT rectangle;
        LPPOINT point = &p;
        LPRECT rect = &rectangle;
        GetCursorPos(point);
        GetWindowRect(hWnd, rect);
        int delta_y = point->x - (rect->right + rect->left) / 2;
        int delta_x= point->y - (rect->bottom + rect->top) / 2;

        DIR = normalize(yRotate(delta_y * angle, DIR));
        DIR = normalize(xRotate(delta_x * angle, DIR));

        SetCursorPos((rect->right + rect->left)/2, (rect->bottom + rect->top)/2);
        InvalidateRect(hWnd, NULL, NULL);
    }
    break;
    case WM_PAINT:
    {
        Draw(hWnd);
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

// Обработчик сообщений для окна "О программе".
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
