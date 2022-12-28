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

#include <qlat-utils/show.h>
#include <qlat-utils/assert.h>

#include <sys/time.h>
#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <set>

#ifdef USE_PAPI
#include <papi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// #if defined USE_MPI || defined USE_QMP
// #define USE_MULTI_NODE
// #endif
// 
// #ifdef USE_MULTI_NODE
// #include <mpi.h>
// #endif

#define TIMER(FNAME)                 \
  static std::string fname = FNAME;  \
  static qlat::Timer timer(fname, true); \
  qlat::TimerCtrl timerctrl(timer);

#define TIMER_VERBOSE(FNAME)         \
  static std::string fname = FNAME;  \
  static qlat::Timer timer(fname, true); \
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

API inline double& get_actual_start_time()
// not affected by Timer::reset()
{
  static double time = get_time();
  return time;
}

inline double get_actual_total_time()
{
  return get_time() - get_actual_start_time();
}

API inline double& get_start_time()
// will be reset by Timer::reset()
{
  static double time = get_actual_start_time();
  return time;
}

inline double get_total_time() { return get_time() - get_start_time(); }

API inline int& get_num_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int num_node = 1;
  return num_node;
}

inline int get_num_node()
{
  return get_num_node_internal();
}

API inline int& get_id_node_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int id_node = 0;
  return id_node;
}

inline int get_id_node()
{
  return get_id_node_internal();
}

inline int get_num_thread()
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

inline int get_id_thread()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

API inline int& get_id_node_in_shuffle_internal()
// initialized in begin_comm in qlat/mpi.h
{
  static int id_node_in_shuffle = 0;
  return id_node_in_shuffle;
}

inline int get_id_node_in_shuffle()
{
  return get_id_node_in_shuffle_internal();
}

inline void display_info(const std::string& str, FILE* fp = NULL)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
    display(str, fp);
  }
}

inline void displayln_info(const std::string& str, FILE* fp = NULL)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
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
  if (0 == get_id_node() && 0 == get_id_thread()) {
    display(minimum_verbose_level, str);
  }
}

inline void displayln_info(const long minimum_verbose_level, const std::string& str)
{
  if (0 == get_id_node() && 0 == get_id_thread()) {
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
                   ssprintf("%s='%s' (default)", var_name.c_str(), x0.c_str()));
    return x0;
  } else {
    displayln_info(0, ssprintf("%s='%s'", var_name.c_str(), val.c_str()));
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
  const long x0 = -1; // default verbose_level
  const std::string var_name = "q_verbose";
  const std::string val = get_env(var_name);
  long x;
  if (val == "") {
    x = x0;
  } else {
    x = read_long(val);
  }
  if (x >= 0) {
    displayln_info(ssprintf("%s=%ld", var_name.c_str(), x));
  }
  return x;
}

API inline long& verbose_level()
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

template <class K>
bool has(const std::set<K>& m, const K& key)
{
  typename std::set<K>::const_iterator it = m.find(key);
  return it != m.end();
}

template <class K, class M>
bool has(const std::map<K, M>& m, const K& key)
{
  typename std::map<K, M>::const_iterator it = m.find(key);
  return it != m.end();
}

API inline void initialize_papi()
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

struct API TimerInfo {
  std::string fname;
  double dtime;
  double accumulated_time;
  long long dflops;
  long long accumulated_flops;
  int call_times;
  //
  TimerInfo() { init(); }
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
  void merge(const TimerInfo& x)
  {
    if (std::isnan(dtime)) {
      dtime = x.dtime;
      dflops = x.dflops;
    }
    accumulated_time += x.accumulated_time;
    accumulated_flops += x.accumulated_flops;
    call_times += x.call_times;
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
    if (0 == get_id_node() && 0 == get_id_thread()) {
      show_avg_always(info, fname_len);
    }
  }
};

inline bool compare_time_info_p(const TimerInfo* p1, const TimerInfo* p2)
{
  return p1->accumulated_time < p2->accumulated_time;
}

struct API Timer {
  const char* cname;
  long info_index;
  bool is_using_total_flops;
  long is_running;
  double start_time;
  double stop_time;
  long long start_flops;
  long long stop_flops;
  long long flops;
  //
  API static std::map<std::string, long>& get_timer_info_index_map()
  {
    static std::map<std::string, long> timer_info_index_map;
    return timer_info_index_map;
  }
  //
  API static std::vector<TimerInfo>& get_timer_database()
  {
    static std::vector<TimerInfo> timer_database;
    return timer_database;
  }
  //
  API static std::vector<std::vector<TimerInfo> >& get_timer_database_history()
  {
    static std::vector<std::vector<TimerInfo> > timer_database_history;
    return timer_database_history;
  }
  //
  API static std::vector<double>& get_start_time_history()
  {
    static std::vector<double> history;
    return history;
  }
  //
  API static std::vector<long>&
  get_max_call_times_for_always_show_info_history()
  {
    static std::vector<long> history;
    return history;
  }
  //
  API static std::vector<long>& get_timer_stack()
  {
    static std::vector<long> stack;
    return stack;
  }
  //
  API static void reset(const long max_call_times_for_always_show_info_ = -1)
  // if max_call_times_for_always_show_info_ <= -1:
  // then do not change the current value.
  // else update the max_call_times_for_always_show_info.
  {
    std::vector<TimerInfo>& tdb = get_timer_database();
    for (long i = 0; i < (long)tdb.size(); ++i) {
      tdb[i].reset();
    }
    get_start_time() = get_time();
    if (max_call_times_for_always_show_info_ >= 0) {
      max_call_times_for_always_show_info() =
          max_call_times_for_always_show_info_;
    }
    std::vector<std::vector<TimerInfo> >& tdb_history =
        get_timer_database_history();
    displayln_info(
        0, ssprintf("Timer::reset(%ld): Reset all timers! (level = %ld)",
                    max_call_times_for_always_show_info(),
                    (long)tdb_history.size()));
  }
  //
  API static void fork(const long max_call_times_for_always_show_info_ = -1)
  // if max_call_times_for_always_show_info_ <= -1:
  // then do not change the current value.
  // else update the max_call_times_for_always_show_info.
  {
    get_start_time_history().push_back(get_start_time());
    get_start_time() = get_time();
    get_max_call_times_for_always_show_info_history().push_back(
        max_call_times_for_always_show_info());
    if (max_call_times_for_always_show_info_ >= 0) {
      max_call_times_for_always_show_info() =
          max_call_times_for_always_show_info_;
    }
    std::vector<std::vector<TimerInfo> >& tdb_history =
        get_timer_database_history();
    std::vector<TimerInfo>& tdb = get_timer_database();
    tdb_history.push_back(tdb);
    for (long i = 0; i < (long)tdb.size(); ++i) {
      tdb[i].reset();
    }
    displayln_info(0,
                   ssprintf("Timer::fork(%ld): Fork all timers! (level = %ld)",
                            max_call_times_for_always_show_info(),
                            (long)tdb_history.size()));
  }
  //
  API static void merge()
  // call merge only after fork
  {
    qassert(get_start_time_history().size() >= 1);
    get_start_time() = get_start_time_history().back();
    get_start_time_history().pop_back();
    qassert(get_max_call_times_for_always_show_info_history().size() >= 1);
    max_call_times_for_always_show_info() =
        get_max_call_times_for_always_show_info_history().back();
    get_max_call_times_for_always_show_info_history().pop_back();
    std::vector<std::vector<TimerInfo> >& tdb_history =
        get_timer_database_history();
    std::vector<TimerInfo>& tdb = get_timer_database();
    qassert(tdb_history.size() >= 1);
    qassert(tdb.size() >= tdb_history.back().size());
    for (long i = 0; i < (long)tdb_history.back().size(); ++i) {
      tdb[i].merge(tdb_history.back()[i]);
    }
    tdb_history.pop_back();
    displayln_info(0,
                   ssprintf("Timer::merge(): Merge all timers! (level = %ld)",
                            (long)tdb_history.size()));
  }
  //
  API static double& minimum_autodisplay_interval()
  // qlat parameter
  {
    static double time =
        get_env_double_default("q_timer_mini_auto_display", 5.0 * 60.0);
    return time;
  }
  //
  API static double& minimum_duration_for_show_info()
  // qlat parameter
  {
    static double time = get_env_double_default("q_timer_mini_auto_show", 1.0);
    return time;
  }
  //
  API static double& minimum_duration_for_show_stop_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  API static double& minimum_duration_for_show_start_info()
  {
    return minimum_duration_for_show_info();
  }
  //
  API static long& max_call_times_for_always_show_info()
  // qlat parameter
  {
    static long max_call_times =
        get_env_long_default("q_timer_max_always_show", 10);
    return max_call_times;
  }
  //
  API static long& max_function_name_length_shown()
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
    is_using_total_flops = false;
    get_start_time();
    initialize_papi();
    info_index = -1;
    is_running = 0;
  }
  void init(const std::string& fname_str)
  {
    std::map<std::string, long>& tiim = get_timer_info_index_map();
    std::vector<TimerInfo>& tdb = get_timer_database();
    if (has(tiim, fname_str)) {
      info_index = tiim[fname_str];
    } else {
      info_index = tdb.size();
      tiim[fname_str] = info_index;
      TimerInfo info;
      info.fname = fname_str;
      tdb.push_back(info);
    }
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
    if (is_running > 0) {
      is_running += 1;
      return;
    } else if (is_running == 0) {
      is_running = 1;
    } else {
      TimerInfo& info = get_timer_database()[info_index];
      info.show_avg_always("debug", max_function_name_length_shown());
      displayln(ssprintf("%s::%s ERROR: is_running=%d", cname,
                         info.fname.c_str(), is_running));
      Timer::display_stack();
      qassert(false);
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
    qassert(not t_stack.empty());
    if (not(t_stack.back() == info_index)) {
      displayln(ssprintf("%s::%s ERROR: stack is corrupted", cname,
                         info.fname.c_str()));
      Timer::display_stack();
      qassert(false);
    }
    t_stack.pop_back();
    if (is_running <= 0) {
      info.show_avg_always("debug", max_function_name_length_shown());
      displayln(ssprintf("%s::%s ERROR: is_running=%d", cname,
                         info.fname.c_str(), is_running));
      Timer::display_stack();
      qassert(false);
    }
    is_running -= 1;
    if (is_running != 0) {
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
          info.dtime >= 2.0 * minimum_duration_for_show_info() ||
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
  API static void test_timer_time_usage()
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
  API static void display(const std::string& tag = "")
  {
    double total_time = get_total_time();
    const std::vector<TimerInfo>& tdb = get_timer_database();
    const std::vector<std::vector<TimerInfo> >& tdb_history =
        get_timer_database_history();
    std::vector<const TimerInfo*> db;
    const long tdbsize = tdb.size();
    for (long i = 0; i < tdbsize; i++) {
      db.push_back(&tdb[i]);
    }
    std::sort(db.begin(), db.end(), compare_time_info_p);
    displayln_info(ssprintf(
        "Timer::display-start: %s (level=%ld) fname : time%% number of calls; "
        "Avg,Tot sec; Avg,Tot flops; Gflops",
        tag.c_str(), (long)tdb_history.size()));
    const long dbsize = db.size();
    for (long i = 0; i < dbsize; i++) {
      if (db[i]->call_times > 0) {
        db[i]->show_avg("display", max_function_name_length_shown());
      }
    }
    displayln_info(ssprintf(
        "Timer::display-end:   %s (level=%ld) --------------------- total %.4E "
        "sec ----------------------",
        tag.c_str(), (long)tdb_history.size(), total_time));
  }
  //
  API static void autodisplay(const double time)
  {
    static double last_time = get_start_time();
    if (time - last_time > minimum_autodisplay_interval()) {
      last_time = time;
      display("autodisplay");
    }
  }
  API static void autodisplay()
  {
    const double time = get_time();
    autodisplay(time);
  }
  //
  API static void display_stack_always()
  // display on any process or thread
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
  API static void display_stack()
  // only display if id_node == 0 and thread_num == 0
  {
    if (0 == get_id_node() && 0 == get_id_thread()) {
      display_stack_always();
    }
  }
};

struct API TimerCtrl {
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
    if (get_id_thread() != 0) return;
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
  qassert(vasprintf(&str, format, args) >= 0);
  display(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

inline void DisplayInfo(const char* cname, const char* fname,
                        const char* format, ...)
{
  int rank = get_id_node();
  if (0 != rank) {
    return;
  }
  va_list args;
  va_start(args, format);
  char* str;
  qassert(vasprintf(&str, format, args) >= 0);
  display_info(ssprintf("%s::%s : %s", cname, fname, str));
  std::free(str);
}

}  // namespace qlat
