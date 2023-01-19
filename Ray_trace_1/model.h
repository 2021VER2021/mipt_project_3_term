#pragma once
#include <ctime>
#include <windows.h>
#include <thread>
#include <atomic>
#include <chrono>

void do_nothing(Render, double) { return; }

class Model {
	time_t last_time;
	time_t cur_time;
	void (&compute)(Render, double);
	std::atomic<bool> RUNNING;
	HWND hWnd = NULL;
public:
	Model() : Model(do_nothing) {}
	Model(void (&f)(Render, double)) : compute(f) {};
	~Model(){}

	void send_wnd(HWND hwnd) {
		hWnd = hwnd;
	}

	void computation(Render r) {
		while (RUNNING.load()) {
			step(r);
			InvalidateRect(hWnd, NULL, NULL);
			std::this_thread::sleep_for(std::chrono::milliseconds(30));
		}
	}

	std::thread start(Render r) {
		time(&last_time);
		RUNNING.store(true);
		return std::thread(&Model::computation, this, r);
	}

	void stop() {
		RUNNING.store(false);
	}

	void step(Render r) {
		time(&cur_time);
		compute(r, difftime(last_time, cur_time));
		last_time = cur_time;
		
		//std::this_thread::sleep_for(std::chrono::milliseconds(100));
		//PostMessageW(hWnd, MODEL_UPDATE, 1, 1);
	}

};
