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

#include <qutils/show.h>

#include <sys/time.h>
#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef USE_PAPI
#include <omp.h>
#include <papi.h>
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

#define TIMER(FNAME)                 \
  static std::string fname = FNAME;  \
  static qlat::Timer timer(fname); \
  qlat::TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE(FNAME)         \
  static std::string fname = FNAME;  \
  static qlat::Timer timer(fname); \
  qlat::TimerCtrl timerctrl(timer, true);

#define TIMER_FLOPS(FNAME)                  \
  static std::string fname = FNAME;         \
  static qlat::Timer timer(fname, false); \
  qlat::TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE_FLOPS(FNAME)        \
  static std::string fname = FNAME;       \
  static qlat::Timer timer(fname, false); \
  qlat::TimerCtrl timerctrl(timer, true);

namespace qlat
{  //

inline long& verbose_level();

inline double get_time()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

inline double& get_actual_start_time()
// not affected by Timer::reset()
{
  static double time = get_time();
  return time;
}

inline double get_actual_total_time()
{
  return get_time() - get_actual_start_time();
}

inline double& get_start_time()
{
  static double time = get_actual_start_time();
  return time;
}

inline double get_total_time() { return get_time() - get_start_time(); }

inline int& get_num_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int num_node = 1;
  return num_node;
}

inline int get_num_node()
{
  return get_num_node_internal();
}

inline int& get_id_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int id_node = 0;
  return id_node;
}

inline int get_id_node()
{
  return get_id_node_internal();
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
  if (0 == get_id_node() && 0 == get_thread_num()) {
    display(str, fp);
  }
}

inline void displayln_info(const std::string& str, FILE* fp = NULL)
{
  if (0 == get_id_node() && 0 == get_thread_num()) {
    displayln(str, fp);
  }
}

inline void display(const long minimum_verbose_level, const std::string& str)
{
  if (verbose_level() >= minimum_verbose_level) {
    display(str);
  }
}

inline void displayln(const long minimum_verbose_level, const std::string& str)
{
  if (verbose_level() >= minimum_verbose_level) {
    displayln(str);
  }
}

inline void display_info(const long minimum_verbose_level, const std::string& str)
{
  if (0 == get_id_node() && 0 == get_thread_num()) {
    display(minimum_verbose_level, str);
  }
}

inline void displayln_info(const long minimum_verbose_level, const std::string& str)
{
  if (0 == get_id_node() && 0 == get_thread_num()) {
    displayln(minimum_verbose_level, str);
  }
}

inline std::string get_env(const std::string& var_name)
{
  const char* value = getenv(var_name.c_str());
  if (value == NULL) {
    return std::string();
  } else {
    return std::string(value);
  }
}

inline std::string get_env_default(const std::string& var_name,
                                   const std::string& x0)
{
  const std::string val = get_env(var_name);
  if (val == "") {
    displayln_info(0,
                   ssprintf("%s=%s (default)", var_name.c_str(), x0.c_str()));
    return x0;
  } else {
    displayln_info(0, ssprintf("%s=%s", var_name.c_str(), val.c_str()));
    return val;
  }
}

inline double get_env_double_default(const std::string& var_name,
                                     const double x0)
{
  const std::string val = get_env(var_name);
  double x;
  if (val == "") {
    x = x0;
    displayln_info(0, ssprintf("%s=%lG (default)", var_name.c_str(), x));
  } else {
    x = read_double(val);
    displayln_info(0, ssprintf("%s=%lG", var_name.c_str(), x));
  }
  return x;
}

inline long get_env_long_default(const std::string& var_name, const long x0)
{
  const std::string val = get_env(var_name);
  long x;
  if (val == "") {
    x = x0;
    displayln_info(0, ssprintf("%s=%ld (default)", var_name.c_str(), x));
  } else {
    x = read_long(val);
    displayln_info(0, ssprintf("%s=%ld", var_name.c_str(), x));
  }
  return x;
}

inline long get_verbose_level()
{
  const long x0 = 0;
  const std::string var_name = "q_verbose";
  const std::string val = get_env(var_name);
  long x;
  if (val == "") {
    x = x0;
    displayln_info(ssprintf("%s=%ld", var_name.c_str(), x));
  } else {
    x = read_long(val);
    if (x >= 0) {
      displayln_info(ssprintf("%s=%ld", var_name.c_str(), x));
    }
  }
  return x;
}

inline long& verbose_level()
// qlat parameter
{
  static long level = get_verbose_level();
  return level;
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
  displayln_info(0, "PAPI::initialize_papi Start.");
  PAPI_library_init(PAPI_VER_CURRENT);
  PAPI_thread_init((unsigned long (*)(void))(omp_get_thread_num));
  initialized = true;
  displayln_info(0, "PAPI::initialize_papi Finish.");
#endif
}

struct TimerInfo {
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
    reset();
  }
  //
  void reset()
  {
    dtime = 0.0 / 0.0;
    accumulated_time = 0;
    dflops = 0;
    accumulated_flops = 0;
    call_times = 0;
  }
  //
  void show_start(const int fname_len) const
  {
    double total_time = get_total_time();
    std::string fnameCut;
    fnameCut.assign(fname, 0, fname_len);
    displayln_info(ssprintf(
        "Timer::start %s :%5.1f%% %8d calls %.3E,%.3E sec %8.3f Gflops (%.3E "
        "flops)",
        ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str())
            .c_str(),
        accumulated_time / total_time * 100, call_times, dtime,
        accumulated_time / (call_times - 1), dflops / dtime / 1.0E9,
        (double)dflops));
  }
  //
  void show_stop(const int fname_len) const
  {
    double total_time = get_total_time();
    std::string fnameCut;
    fnameCut.assign(fname, 0, fname_len);
    displayln_info(ssprintf(
        "Timer::stop  %s :%5.1f%% %8d calls %.3E,%.3E sec %8.3f Gflops (%.3E "
        "flops)",
        ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str())
            .c_str(),
        accumulated_time / total_time * 100, call_times, dtime,
        accumulated_time / call_times, dflops / dtime / 1.0E9, (double)dflops));
  }
  //
  void show_avg_always(const std::string& info, const int fname_len) const
  {
    double total_time = get_total_time();
    std::string fnameCut;
    fnameCut.assign(fname, 0, fname_len);
    displayln(ssprintf(
        "Timer::%s %s :%7.3f%% %8d calls; %.2E,%.2E sec; %.2E,%.2E flops; "
        "%5.2f Gflops",
        info.c_str(),
        ssprintf(ssprintf("%%%ds", fname_len).c_str(), fnameCut.c_str())
            .c_str(),
        accumulated_time / total_time * 100, call_times,
        accumulated_time / call_times, accumulated_time,
        (double)accumulated_flops / (double)call_times,
        (double)accumulated_flops,
        accumulated_flops / accumulated_time / 1.0E9));
  }
  void show_avg(const std::string& info, const int fname_len) const
  {
    if (0 == get_id_node() && 0 == get_thread_num()) {
      show_avg_always(info, fname_len);
    }
  }
};

inline bool compare_time_info_p(const TimerInfo* p1, const TimerInfo* p2)
{
  return p1->accumulated_time < p2->accumulated_time;
}

struct Timer {
  const char* cname;
  long info_index;
  bool is_using_total_flops;
  long isRunning;
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
  static std::vector<long>& get_timer_stack()
  {
    static std::vector<long> stack;
    return stack;
  }
  //
  static void reset()
  {
    displayln_info(0, "Timer::reset(): Reset all timers!");
    std::vector<TimerInfo>& tdb = get_timer_database();
    for (long i = 0; i < (long)tdb.size(); ++i) {
      tdb[i].reset();
    }
    get_start_time() = get_time();
  }
  //
  static double& minimum_autodisplay_interval()
  // qlat parameter
  {
    static double time =
        get_env_double_default("q_timer_mini_auto_display", 5.0 * 60.0);
    return time;
  }
  //
  static double& minimum_duration_for_show_info()
  // qlat parameter
  {
    static double time = get_env_double_default("q_timer_mini_auto_show", 1.0);
    return time;
  }
  //
  static double& minimum_duration_for_show_stop_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  static double& minimum_duration_for_show_start_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  static long& max_call_times_for_always_show_info()
  // qlat parameter
  {
    static long max_call_times =
        get_env_long_default("q_timer_max_always_show", 10);
    return max_call_times;
  }
  //
  static long& max_function_name_length_shown()
  // qlat parameter
  {
    static long max_len = get_env_long_default("q_timer_max_func_name_len", 50);
    return max_len;
  }
  //
  Timer() { init(); }
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
    isRunning = 0;
  }
  void init(const std::string& fname_str)
  {
    std::vector<TimerInfo>& tdb = get_timer_database();
    const long size = tdb.size();
    for (long i = 0; i < size; i++) {
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
    get_timer_stack().push_back(info_index);
    if (isRunning > 0) {
      isRunning += 1;
      return;
    } else if (isRunning == 0) {
      isRunning = 1;
    } else {
      TimerInfo& info = get_timer_database()[info_index];
      info.show_avg_always("debug", max_function_name_length_shown());
      displayln(ssprintf("%s::%s ERROR: isRunning=%d", cname,
                         info.fname.c_str(), isRunning));
      Timer::display_stack();
      pqassert(false);
    }
    TimerInfo& info = get_timer_database()[info_index];
    info.call_times++;
    if (verbose_level() > 0) {
      if (verbose ||
          info.accumulated_time >=
              info.call_times * minimum_duration_for_show_info() ||
          info.call_times <= max_call_times_for_always_show_info()) {
        info.show_start(max_function_name_length_shown());
      }
    }
    start_flops = is_using_total_flops ? get_total_flops() : 0;
    flops = 0;
    start_time = get_time();
  }
  //
  void stop(bool verbose = false)
  {
    TimerInfo& info = get_timer_database()[info_index];
    std::vector<long>& t_stack = get_timer_stack();
    pqassert(not t_stack.empty());
    if (not (t_stack.back() == info_index)) {
      displayln(ssprintf("%s::%s ERROR: stack is corrupted", cname,
                         info.fname.c_str()));
      Timer::display_stack();
      pqassert(false);
    }
    t_stack.pop_back();
    if (isRunning <= 0) {
      info.show_avg_always("debug", max_function_name_length_shown());
      displayln(ssprintf("%s::%s ERROR: isRunning=%d", cname,
                         info.fname.c_str(), isRunning));
      Timer::display_stack();
      pqassert(false);
    }
    isRunning -= 1;
    if (isRunning != 0) {
      return;
    }
    stop_time = get_time();
    if (is_using_total_flops) {
      stop_flops = get_total_flops();
    } else {
      stop_flops = start_flops + flops;
    }
    info.dtime = stop_time - start_time;
    info.dflops = stop_flops - start_flops;
    bool is_show = false;
    if (verbose_level() > 0) {
      if (verbose ||
          info.accumulated_time >=
              info.call_times * minimum_duration_for_show_info() ||
          info.dtime >= 5.0 * minimum_duration_for_show_info() ||
          info.call_times <= max_call_times_for_always_show_info()) {
        is_show = true;
      }
    }
    info.accumulated_time += info.dtime;
    info.accumulated_flops += info.dflops;
    if (is_show) {
      info.show_stop(max_function_name_length_shown());
    }
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
    const long tdbsize = tdb.size();
    for (long i = 0; i < tdbsize; i++) {
      db.push_back(&tdb[i]);
    }
    std::sort(db.begin(), db.end(), compare_time_info_p);
    displayln_info(
        ssprintf("Timer::display-start: %s fname : time%% number of calls; "
                 "Avg,Tot sec; Avg,Tot flops; Gflops",
                 str.c_str()));
    const long dbsize = db.size();
    for (long i = 0; i < dbsize; i++) {
      db[i]->show_avg("display", max_function_name_length_shown());
    }
    displayln_info(
        ssprintf("Timer::display-end:   %s --------------------- total %.4E "
                 "sec ----------------------",
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
  //
  static void display_stack_always()
  {
    displayln("display_stack start");
    const std::vector<TimerInfo>& tdb = get_timer_database();
    const std::vector<long>& t_stack = get_timer_stack();
    for (long i = (long)t_stack.size() - 1; i >= 0; --i) {
      const long info_index = t_stack[i];
      tdb[info_index].show_avg_always(ssprintf("stack[%3ld]", i),
                                      max_function_name_length_shown());
    }
    displayln("display_stack end");
  }
  static void display_stack()
  {
    if (0 == get_id_node() && 0 == get_thread_num()) {
      display_stack_always();
    }
  }
};

struct TimerCtrl {
  Timer* ptimer;
  bool verbose;
  //
  TimerCtrl() { init(); }
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

inline void* tmalloc(size_t size) { return timer_malloc(size); }

inline void tfree(void* ptr) { timer_free(ptr); }

///////////////////////////////////////////////////////////////////////

inline void Display(const char* cname, const char* fname, const char* format,
                    ...)
{
  va_list args;
  va_start(args, format);
  char* str;
  pqassert(vasprintf(&str, format, args) >= 0);
  display(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

inline void DisplayInfo(const char* cname, const char* fname,
                        const char* format, ...)
{
  static int rank = get_id_node();
  if (0 != rank) {
    return;
  }
  va_list args;
  va_start(args, format);
  char* str;
  pqassert(vasprintf(&str, format, args) >= 0);
  display_info(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

}  // namespace qlat

#ifndef USE_NAMESPACE
using namespace qlat;
#endif
