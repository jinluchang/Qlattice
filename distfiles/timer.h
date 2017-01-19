// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2014 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include "show.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>

#ifdef USE_PAPI
#include <papi.h>
#include <omp.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined USE_MPI || defined USE_QMP
#define USE_MULTI_NODE
#endif

#ifdef USE_MULTI_NODE
#include <mpi.h>
#endif

#define TIMER(FNAME) \
  static const char* fname = FNAME; \
  static qtimer::Timer timer(fname); \
  qtimer::TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE(FNAME) \
  static const char* fname = FNAME; \
  static qtimer::Timer timer(fname); \
  qtimer::TimerCtrl timerctrl(timer, true);

#define TIMER_FLOPS(FNAME) \
  static const char* fname = FNAME; \
  static qtimer::Timer timer(fname, false); \
  qtimer::TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE_FLOPS(FNAME) \
  static const char* fname = FNAME; \
  static qtimer::Timer timer(fname, false); \
  qtimer::TimerCtrl timerctrl(timer, true);

namespace qtimer {

using namespace qshow;

inline double get_time()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

inline double get_start_time()
{
  static double time = get_time();
  return time;
}

inline double get_total_time()
{
  return get_time() - get_start_time();
}

inline int compute_rank()
{
#ifdef USE_MULTI_NODE
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  return myid;
#else
  return 0;
#endif
}

inline int get_rank()
{
  static int myid = compute_rank();
  return myid;
}

inline int get_thread_num()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline void display_info(const std::string& str, FILE* fp = NULL)
{
  if (0 == get_rank() && 0 == get_thread_num()) {
    display(str, fp);
  }
}

inline void displayln_info(const std::string& str, FILE* fp = NULL)
{
  if (0 == get_rank() && 0 == get_thread_num()) {
    displayln(str, fp);
  }
}

inline long long get_total_flops()
{
  long long flops = 0;
#ifdef USE_PAPI
  const int n_threads = omp_get_max_threads();
  long long flopses[n_threads];
  std::memset(flopses, 0, n_threads * sizeof(long long));
#pragma omp parallel
  {
    float rtime, ptime, mflops;
    int i = omp_get_thread_num();
    PAPI_flops(&rtime, &ptime, &flopses[i], &mflops);
  }
  for (int i = 0; i < n_threads; i++) {
    flops += flopses[i];
  }
#endif
  return flops;
}

inline void initialize_papi()
{
#ifdef USE_PAPI
  static bool initialized = false;
  if (initialized) {
    return;
  }
  displayln_info("PAPI::initialize_papi Start.");
  PAPI_library_init(PAPI_VER_CURRENT);
  PAPI_thread_init((unsigned long (*)(void)) (omp_get_thread_num));
  initialized = true;
  displayln_info("PAPI::initialize_papi Finish.");
#endif
}

struct TimerInfo
{
  std::string fname;
  double dtime;
  double accumulated_time;
  long long dflops;
  long long accumulated_flops;
  int call_times;
  //
  void init()
  {
    fname = "Unknown";
    dtime = 0.0 / 0.0;
    accumulated_time = 0;
    dflops = 0;
    accumulated_flops = 0;
    call_times = 0;
  }
  //
  void show_last(const char* info, const int fname_len) const
  {
    double total_time = get_total_time();
    std::string fnameCut;
    fnameCut.assign(fname, 0, fname_len);
    displayln_info(
        ssprintf("Timer::%s %s :%5.1f%% %8d calls %.3E sec %8.3f Gflops (%.3E flops)",
          NULL == info ? "" : info,
          ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str()).c_str(),
          accumulated_time / total_time * 100, call_times,
          dtime,
          dflops / dtime / 1.0E9,
          (double)dflops));
  }
  //
  void show_avg(const char* info, const int fname_len) const
  {
    double total_time = get_total_time();
    std::string fnameCut;
    fnameCut.assign(fname, 0, fname_len);
    displayln_info(
        ssprintf("Timer::%s %s :%7.3f%% %8d calls; %.2E,%.2E sec; %.2E,%.2E flops; %5.2f Gflops",
          NULL == info ? "" : info,
          ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str()).c_str(),
          accumulated_time / total_time * 100, call_times,
          accumulated_time / call_times,
          accumulated_time,
          (double)accumulated_flops / (double)call_times,
          (double)accumulated_flops,
          accumulated_flops / accumulated_time / 1.0E9));
  }
};

inline bool compare_time_info_p(const TimerInfo* p1, const TimerInfo* p2)
{
  return p1->accumulated_time < p2->accumulated_time;
}

struct Timer
{
  const char* cname;
  int info_index;
  bool is_using_total_flops;
  bool isRunning;
  double start_time;
  double stop_time;
  long long start_flops;
  long long stop_flops;
  long long flops;
  //
  static std::vector<TimerInfo>& get_timer_database()
  {
    static std::vector<TimerInfo> timer_database;
    return timer_database;
  }
  //
  static double& minimum_autodisplay_interval()
  {
    static double time = 365 * 24 * 3600.0;
    return time;
  }
  //
  static double& minimum_duration_for_show_stop_info()
  {
    static double time = 60.0;
    return time;
  }
  //
  static double& minimum_duration_for_show_start_info()
  {
    static double time = 60.0;
    return time;
  }
  //
  static int& max_call_times_for_always_show_info()
  {
    static int max_call_times = 10;
    return max_call_times;
  }
  //
  static int& max_function_name_length_shown()
  {
    static int max_len = 30;
    return max_len;
  }
  //
  Timer()
  {
    init();
  }
  Timer(const std::string& fname_str)
  {
    init();
    init(fname_str);
  }
  Timer(const std::string& cname_str, const std::string& fname_str)
  {
    init();
    init(cname_str, fname_str);
  }
  Timer(const std::string& fname_str, const bool is_using_total_flops_)
  {
    init();
    init(fname_str);
    is_using_total_flops = is_using_total_flops_;
  }
  //
  void init()
  {
    cname = "Timer";
    is_using_total_flops = true;
    get_start_time();
    initialize_papi();
    info_index = -1;
    isRunning = false;
  }
  void init(const std::string& fname_str)
  {
    std::vector<TimerInfo>& tdb = get_timer_database();
    const int size = tdb.size();
    for (int i = 0; i < size; i++) {
      if (fname_str == tdb[i].fname) {
        info_index = i;
        return;
      }
    }
    info_index = tdb.size();
    TimerInfo info;
    tdb.push_back(info);
    tdb[info_index].init();
    tdb[info_index].fname = fname_str;
  }
  void init(const std::string& cname_str, const std::string& fname_str)
  {
    std::string fname = "";
    fname += cname_str;
    fname += "::";
    fname += fname_str;
    init(fname);
  }
  //
  void start(bool verbose = false)
  {
    if (isRunning) {
      return;
    } else {
      isRunning = true;
    }
    TimerInfo& info = get_timer_database()[info_index];
    info.call_times++;
    if (verbose || info.call_times <= max_call_times_for_always_show_info() || info.dtime >= minimum_duration_for_show_start_info()) {
      info.show_last("start", max_function_name_length_shown());
    }
    start_flops = is_using_total_flops ? get_total_flops() : 0 ;
    flops = 0;
    start_time = get_time();
  }
  //
  void stop(bool verbose = false)
  {
    stop_time = get_time();
    assert(isRunning);
    isRunning = false;
    if (is_using_total_flops) {
      stop_flops = get_total_flops();
    } else {
      stop_flops = start_flops + flops;
    }
    TimerInfo& info = get_timer_database()[info_index];
    info.dtime = stop_time - start_time;
    info.dflops = stop_flops - start_flops;
    info.accumulated_time += info.dtime;
    info.accumulated_flops += info.dflops;
    if (verbose || info.call_times <= max_call_times_for_always_show_info() || info.dtime >= minimum_duration_for_show_stop_info()) {
      info.show_last("stop ", max_function_name_length_shown());
    }
    autodisplay(stop_time);
  }
  //
  static void test_timer_time_usage()
  {
    {
      static Timer timer("Timer", false);
      static Timer timer_test("Timer-test");
      timer.start();
      timer_test.start();
      timer_test.stop();
      timer.stop();
    }
    {
      static Timer timer_noflop("Timer-noflop", false);
      static Timer timer_test("Timer-test", false);
      timer_noflop.start();
      timer_test.start();
      timer_test.stop();
      timer_noflop.stop();
    }
  }
  //
  static void display(const std::string& str = "")
  {
    double total_time = get_total_time();
    const std::vector<TimerInfo>& tdb = get_timer_database();
    std::vector<const TimerInfo*> db;
    const int tdbsize = tdb.size();
    for (int i = 0; i < tdbsize; i++) {
      db.push_back(&tdb[i]);
    }
    std::sort(db.begin(), db.end(), compare_time_info_p);
    displayln_info(ssprintf(
          "Timer::display-start: %s fname : time%% number of calls; Avg,Tot sec; Avg,Tot flops; Gflops",
          str.c_str()));
    const int dbsize = db.size();
    for (int i = 0; i < dbsize; i++) {
      db[i]->show_avg("display", max_function_name_length_shown());
    }
    displayln_info(ssprintf(
          "Timer::display-end:   %s --------------------- total %.4E sec ----------------------",
          str.c_str(), total_time));
  }
  //
  static void autodisplay(const double time)
  {
    static double last_time = get_start_time();
    if (time - last_time > minimum_autodisplay_interval()) {
      last_time = time;
      display("autodisplay");
    }
  }
  static void autodisplay()
  {
    const double time = get_time();
    autodisplay(time);
  }
};

struct TimerCtrl
{
  Timer* ptimer;
  bool verbose;
  //
  TimerCtrl()
  {
    init();
  }
  TimerCtrl(Timer& timer, bool verbose_ = false)
  {
    init();
    init(timer, verbose_);
  }
  //
  ~TimerCtrl()
  {
    if (NULL != ptimer) {
      ptimer->stop(verbose);
    }
  }
  //
  void init()
  {
    ptimer = NULL;
    verbose = false;
  }
  void init(Timer& timer, bool verbose_ = false)
  {
    if (get_thread_num() != 0) return;
    ptimer = &timer;
    verbose = verbose_;
    ptimer->start(verbose);
  }
};

inline void* timer_malloc(size_t size)
{
  TIMER_FLOPS("timer_malloc");
  timer.flops += size;
  void* p = malloc(size);
  memset(p, 0, size);
  return p;
}

inline void timer_free(void* ptr)
{
  TIMER_FLOPS("timer_free");
  free(ptr);
}

inline void* tmalloc(size_t size)
{
  return timer_malloc(size);
}

inline void tfree(void* ptr)
{
  timer_free(ptr);
}

///////////////////////////////////////////////////////////////////////

inline void Display(const char* cname, const char* fname, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  char* str;
  vasprintf(&str, format, args);
  display(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

inline void DisplayInfo(const char* cname, const char* fname, const char* format, ...)
{
  static int rank = get_rank();
  if (0 != rank) {
    return;
  }
  va_list args;
  va_start(args, format);
  char* str;
  vasprintf(&str, format, args);
  display_info(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

}

#ifndef USE_NAMESPACE
using namespace qtimer;
#endif
