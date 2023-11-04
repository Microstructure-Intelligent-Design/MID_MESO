#pragma once
#include "time.h"
#include<iostream>
#include<string>
namespace pf {
	using namespace std;
	struct ints_time {
		int year;
		int month;
		int day;
		int hour;
		int minute;
		int second;
		ints_time() {
			year = 0;
			month = 0;
			day = 0;
			hour = 0;
			minute = 0;
			second = 0;
		}
		ints_time& operator=(const ints_time& n) {
			year = n.year;
			month = n.month;
			day = n.day;
			hour = n.hour;
			minute = n.minute;
			second = n.second;
			return *this;
		}
		ints_time& operator+(const ints_time& n) {
			year += n.year;
			month += n.month;
			day += n.day;
			hour += n.hour;
			minute += n.minute;
			second += n.second;
			return *this;
		}
		ints_time& operator-(const ints_time& n) {
			year -= n.year;
			month -= n.month;
			day -= n.day;
			hour -= n.hour;
			minute -= n.minute;
			second -= n.second;
			return *this;
		}
		bool operator>(const ints_time& n) {
			if (year > n.year)
				return true;
			else if (year < n.year)
				return false;
			if (month > n.month)
				return true;
			else if (month < n.month)
				return false;
			if (day > n.day)
				return true;
			else if (day < n.day)
				return false;
			if (hour > n.hour)
				return true;
			else if (hour < n.hour)
				return false;
			if (minute > n.minute)
				return true;
			else if (minute < n.minute)
				return false;
			if (second > n.second)
				return true;
			else if (second < n.second)
				return false;
			return false;
		}
		bool operator==(const ints_time& n) {
			if (year == n.year && month == n.month && day == n.day && hour == n.hour && minute == n.minute && second == n.second)
				return true;
			return false;
		}
	};
	namespace timer{
		static void init();
		static double total_duration_sec();
		static double total_duration_min();
		static double total_duration_hour();
		static void print_cunrrent_time_on_screen(); 
		static string return_cunrrent_time_by_string();
		static string get_cunrrent_time_as_file_name();
		static void get_cunrrent_time_by_int(int& year, int& month, int& day, int& hour, int& minute, int& second);
		static void get_cunrrent_time_by_ints_time(ints_time& t);
		static void lock_class();
		static void unlock_class();
		static void interval_begin();
		static double interval_end();  //By Second

		static bool classlock = false;
		static clock_t init_t;
		static clock_t interval_t;
	};
	inline double timer::interval_end() {
		return (double)((clock() - interval_t) / CLOCKS_PER_SEC);
	}
	inline void timer::interval_begin() {
		interval_t = clock();
	}
	inline void timer::init(){
		init_t = clock();
	}
	inline double timer::total_duration_sec(){
		return (double)((clock() - init_t) / CLOCKS_PER_SEC);
	}
	inline double timer::total_duration_min() {
		return (double)((clock() - init_t) / CLOCKS_PER_SEC / 60.0);
	}
	inline double timer::total_duration_hour() {
		return (double)((clock() - init_t) / CLOCKS_PER_SEC / 3600.0);
	}
	inline void timer::lock_class() {
		classlock = true;
	}
	inline void timer::unlock_class() {
		classlock = false;
	}
	inline void timer::print_cunrrent_time_on_screen() {
		if (classlock == false) {
			time_t timer;
			time(&timer);
			#ifdef _WIN32
				char buf[26];
				ctime_s(buf, 26, &timer);                                     ///< Windows style directory separator
			#else
				char* buf;
				buf = ctime(&timer);                                         ///< Unix/Linux style directory separator
			#endif
			std::cout << "## " << buf;
		}
	}
	inline string timer::return_cunrrent_time_by_string() {
		if (classlock == false) {
			time_t timer;
			time(&timer);
#ifdef _WIN32
			char buf[26];
			ctime_s(buf, 26, &timer);                                     ///< Windows style directory separator
#else
			char* buf;
			buf = ctime(&timer);                                         ///< Unix/Linux style directory separator
#endif
			string time(buf);
			return "## " + time;
		}
		return "";
	}
	inline string timer::get_cunrrent_time_as_file_name() {
		if (classlock == false) {
			time_t timer;
			time(&timer);
#ifdef _WIN32
			char buf[26];
			ctime_s(buf, 26, &timer);                                     ///< Windows style directory separator
#else
			char* buf;
			buf = ctime(&timer);                                         ///< Unix/Linux style directory separator
#endif
			string time(buf);

			for (auto charr = time.begin(); charr < time.end();) {
				if (*charr == '\n')
					time.erase(charr);
				else if (*charr == ':') {
					*charr = '-';
					charr++;
				}
				else
					charr++;
			}

			return time;
		}
		return "";
	}
	inline void timer::get_cunrrent_time_by_int(int& year, int& month, int& day, int& hour, int& minute, int& second) {
		const time_t t = time(NULL);
		struct tm* system_time = localtime(&t);
		year = 1900 + system_time->tm_year;
		month = 1 + system_time->tm_mon;
		day = system_time->tm_mday;
		hour = system_time->tm_hour;
		minute = system_time->tm_min;
		second = system_time->tm_sec;
	}
	inline void timer::get_cunrrent_time_by_ints_time(ints_time& tt) {
		const time_t t = time(NULL);
		struct tm* system_time = localtime(&t);
		tt.year = 1900 + system_time->tm_year;
		tt.month = 1 + system_time->tm_mon;
		tt.day = system_time->tm_mday;
		tt.hour = system_time->tm_hour;
		tt.minute = system_time->tm_min;
		tt.second = system_time->tm_sec;
	}
}