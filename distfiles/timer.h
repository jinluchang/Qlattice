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

#ifndef INCLUDED_TIMER_H
#define INCLUDED_TIMER_H

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

#if defined USE_MPI || defined USE_QMP
#define USE_MULTI_NODE
#endif

#ifdef USE_MULTI_NODE
#include <mpi.h>
#endif

inline double getTime()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

inline int getRank()
{
#ifdef USE_MULTI_NODE
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  return myid;
#else
  return 0;
#endif
}

inline double getStartTime()
{
  static double time = getTime();
  return time;
}

inline double getTotalTime()
{
  return getTime() - getStartTime();
}

inline void Display(const char* cname, const char* fname, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  char* str;
  vasprintf(&str, format, args);
  std::printf("%s::%s : %s", cname, fname, str);
  std::free(str);
}

inline void DisplayInfo(const char* cname, const char* fname, const char* format, ...)
{
  static int rank = getRank();
  if (0 != rank) {
    return;
  }
  va_list args;
  va_start(args, format);
  char* str;
  vasprintf(&str, format, args);
  std::printf("%s::%s : %s", cname, fname, str);
  std::free(str);
}

inline long long getTotalFlops()
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

inline void initializePAPI()
{
#ifdef USE_PAPI
  static bool initialized = false;
  if (initialized) {
    return;
  }
  DisplayInfo("PAPI", "initializePAPI", "Start.\n");
  PAPI_library_init(PAPI_VER_CURRENT);
  PAPI_thread_init((unsigned long (*)(void)) (omp_get_thread_num));
  initialized = true;
  DisplayInfo("PAPI", "initializePAPI", "Finish.\n");
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
  void showLast(const char* info = NULL) const
  {
    double total_time = getTotalTime();
    std::string fnameCut;
    fnameCut.assign(fname, 0, 30);
    DisplayInfo("Timer", NULL == info ? "" : info,
        "%30s :%5.1f%%%9d calls. Last %.3E secs%8.3f Gflops (%.3E flops)\n",
        fnameCut.c_str(),
        accumulated_time / total_time * 100, call_times,
        dtime,
        dflops / dtime / 1.0E9,
        (double)dflops);
  }
  //
  void showAvg(const char* info = NULL) const
  {
    double total_time = getTotalTime();
    std::string fnameCut;
    fnameCut.assign(fname, 0, 30);
    DisplayInfo("Timer", NULL == info ? "" : info,
        "%30s :%7.3f%%%9d calls; %.2E,%.2E secs; %.2E,%.2E flops;%6.2f Gflops\n",
        fnameCut.c_str(),
        accumulated_time / total_time * 100, call_times,
        accumulated_time / call_times,
        accumulated_time,
        (double)accumulated_flops / (double)call_times,
        (double)accumulated_flops,
        accumulated_flops / accumulated_time / 1.0E9);
  }
  //
  void show(const char* info = NULL) const
  {
    showAvg(info);
  }
};

inline bool compareTimeInfoP(const TimerInfo* p1, const TimerInfo* p2)
{
  return p1->accumulated_time < p2->accumulated_time;
}

struct Timer {
  const char* cname;
  int info_index;
  bool isUsingTotalFlops;
  bool isRunning;
  double start_time;
  double stop_time;
  long long start_flops;
  long long stop_flops;
  long long flops;
  //
  static std::vector<TimerInfo>& getTimerDatabase()
  {
    static std::vector<TimerInfo> timerDatabase;
    return timerDatabase;
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
  Timer(const std::string& fname_str, const bool isUsingTotalFlops_)
  {
    init();
    init(fname_str);
    isUsingTotalFlops = isUsingTotalFlops_;
  }
  //
  void init()
  {
    cname = "Timer";
    isUsingTotalFlops = true;
    getStartTime();
    initializePAPI();
    info_index = -1;
    isRunning = false;
  }
  void init(const std::string& fname_str)
  {
    std::vector<TimerInfo>& tdb = getTimerDatabase();
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
    TimerInfo& info = getTimerDatabase()[info_index];
    info.call_times++;
    if (verbose || info.call_times == 1 || info.dtime >= minimum_duration_for_show_start_info()) {
      info.showLast("start");
    }
    start_flops = isUsingTotalFlops ? getTotalFlops() : 0 ;
    flops = 0;
    start_time = getTime();
  }
  //
  void stop(bool verbose = false)
  {
    stop_time = getTime();
    assert(isRunning);
    isRunning = false;
    if (isUsingTotalFlops) {
      stop_flops = getTotalFlops();
    } else {
      stop_flops = start_flops + flops;
    }
    TimerInfo& info = getTimerDatabase()[info_index];
    info.dtime = stop_time - start_time;
    info.dflops = stop_flops - start_flops;
    info.accumulated_time += info.dtime;
    info.accumulated_flops += info.dflops;
    if (verbose || info.call_times == 1 || info.dtime >= minimum_duration_for_show_stop_info()) {
      info.showLast("stop ");
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
    double total_time = getTotalTime();
    const std::vector<TimerInfo>& tdb = getTimerDatabase();
    std::vector<const TimerInfo*> db;
    const int tdbsize = tdb.size();
    for (int i = 0; i < tdbsize; i++) {
      db.push_back(&tdb[i]);
    }
    std::sort(db.begin(), db.end(), compareTimeInfoP);
    DisplayInfo("Timer", "display-start", "%s fname : time%% number of calls; Avg,Tot secs; Avg,Tot flops; Gflops\n", str.c_str(), total_time);
    const int dbsize = db.size();
    for (int i = 0; i < dbsize; i++) {
      db[i]->showAvg("display");
    }
    DisplayInfo("Timer", "display-end  ", "%s --------------------- total %.4E secs ----------------------\n", str.c_str(), total_time);
  }
  //
  static void autodisplay(const double time)
  {
    static double last_time = getStartTime();
    if (time - last_time > minimum_autodisplay_interval()) {
      last_time = time;
      display("autodisplay");
    }
  }
  static void autodisplay()
  {
    const double time = getTime();
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
    ptimer = NULL;
    verbose = false;
  }
  TimerCtrl(Timer& timer, bool verbose_ = false)
  {
    init(timer, verbose_);
  }
  //
  ~TimerCtrl()
  {
    ptimer->stop(verbose);
  }
  //
  void init(Timer& timer, bool verbose_ = false)
  {
    ptimer = &timer;
    verbose = verbose_;
    ptimer->start(verbose);
  }
};

#define TIMER(FNAME) \
  static const char* fname = FNAME; \
  static Timer timer(fname); \
  TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE(FNAME) \
  static const char* fname = FNAME; \
  static Timer timer(fname); \
  TimerCtrl timerctrl(timer, true);

#define TIMER_FLOPS(FNAME) \
  static const char* fname = FNAME; \
  static Timer timer(fname, false); \
  TimerCtrl timerctrl(timer); \

#define TIMER_VERBOSE_FLOPS(FNAME) \
  static const char* fname = FNAME; \
  static Timer timer(fname, false); \
  TimerCtrl timerctrl(timer, true); \

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
  TIMER("timer_free");
  free(ptr);
}

#define tmalloc(x) timer_malloc(x)

#define tfree(x) timer_free(x)

#endif
